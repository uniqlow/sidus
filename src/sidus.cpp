/**
 * sidus - a converter for Yale Bright Star type catalogs
 *
 * Copyright (c) 2019 Jon Olsson <jlo@wintermute.net>
 *
 * cf. http://tdc-www.harvard.edu/catalogs/catalogsb.html
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#if defined(__unix__)
#include <sys/stat.h>
#include <fcntl.h>
#endif
#include <string>
#include <map>
#include <fstream>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cfloat>
#include <cassert>

namespace {

enum class Epoch { AUTO, J2000, B1950 };
enum class Endian { AUTO, LITTLE, BIG };

struct Header {
	int numStars;
	enum StarId {
		NO_STAR_ID,
		CATALOG_STAR_ID,
		GSC_STAR_ID,
		TYCHO_STAR_ID,
		INTEGER_STAR_ID
	} starId;
	int starNameLength;
	enum ProperMotion {
		NO_PROPER_MOTION,
		PROPER_MOTION,
		RADIAL_VELOCITY
	} properMotion;
	int numMagnitudes;
	int apparentMagnitude;
	int numBytesPerStar;
	Epoch epoch;
	bool littleEndian;
};

struct Star {
	std::string name;
	double rightAscension;	// J2000 or B1950, radians
	double declination;	// J2000 or B1950, radians
	double starId;
	float magnitude;
	struct ProperMotion {
		float rightAscension;	// Radians per year
		float declination;	// Radians per year
	} properMotion;
	double radialVelocity;		// kilometers per second
	char spectralType[3];
};

static
void
usage(FILE * f)
{
	std::fprintf(f, "Usage: sidus [option(s)] <input-file>\n");
	std::fprintf(f, "Options:\n");
	std::fprintf(f, " -a<0-9>	specify apparent magnitude, if multiple exist\n");
	std::fprintf(f, " -f<0-9>	filter magnitudes weaker than specified\n");
	std::fprintf(f, " -B1950		expect B1950 epoch\n");
	std::fprintf(f, " -J2000		expect J2000 epoch\n");
	std::fprintf(f, " -c		output a C header instead of a CSV text\n");
	std::fprintf(f, " -le		expect little-endian format (default)\n");
	std::fprintf(f, " -be		expect big-endian format\n");
	std::fprintf(f, " -s		output single-precision floating point\n");
	std::fprintf(f, " -i		output only information from catalog header\n");
	std::fprintf(f, " -m		sort output by decreasing magnitude\n");
	std::fprintf(f, " -r		sort output by increasing right-ascension\n");
	std::fprintf(f, " -n		output star names\n");
	std::fprintf(f, " -p		output spectral class\n");
	std::fprintf(f, " -h | --help	show this help information\n");
	std::fprintf(f, " -v | --version	show version information\n");
}

static
void
version()
{
	std::fprintf(stdout, "sidus v0.1 by Jon Olsson <jlo@wintermute.net>\n");
}

static
int
getFileSize(size_t* size, const char* const path)
{
#if defined(__unix__)
	struct stat sb;
	if (stat(path, &sb) != -1) {
		*size = sb.st_size;
		return 0;
	}
	return -1;
#else
	FILE *f = fopen(path, "rb");
	if (f) {
		if (fseek(f, 0, SEEK_END) != -1) {
			long rv = ftell(f);
			if (rv >= 0) {
				*size = (size_t)rv;
				fclose(f);
				return 0;
			}
		}
		fclose(f);
	}
	return -1;
#endif
}

static
int
readFile(void * data, char const * const path, size_t const filesize)
{
	unsigned char * buf = (unsigned char*)data;
	FILE *f = fopen(path, "rb");
	if (f) {
		for (size_t i = 0u; i < filesize;) {
			size_t rv = fread(buf + i, 1, filesize - i, f);
			if (rv > 0) {
				i += rv;
			}
			else if (ferror(f)) {
				fclose(f);
				return -2;
			}
		}
		fclose(f);
		return 0;
	}
	return -1;
}


static
void
parse(std::int16_t* val, unsigned char const* const data, bool const littleEndian)
{
	if (littleEndian) {
		*val =
		    ((short)(data[0])<<0) |
		    ((short)(data[1])<<8);
	} else {
		*val = 
		    ((short)(data[0])<<8) |
		    ((short)(data[1])<<0);
	}
}

static
void
parse(std::int32_t* val, unsigned char const* const data, bool const littleEndian)
{
	if (littleEndian) {
		*val =
		    ((int)(data[0])<<0) |
		    ((int)(data[1])<<8) |
		    ((int)(data[2])<<16) |
		    ((int)(data[3])<<24);
	} else {
		*val =
		    ((int)(data[0])<<24) |
		    ((int)(data[1])<<16) |
		    ((int)(data[2])<<8) |
		    ((int)(data[3])<<0);
	}
}

static
void
parse(float* val, unsigned char const* const data, bool const littleEndian)
{
	assert(sizeof(float) == sizeof(std::uint32_t));
	if (littleEndian) {
		auto const v2 =
		    ((std::uint32_t)(data[0])<<0) |
		    ((std::uint32_t)(data[1])<<8) |
		    ((std::uint32_t)(data[2])<<16) |
		    ((std::uint32_t)(data[3])<<24);
		std::memcpy(val, &v2, sizeof v2);
	} else {
		auto const v2 =
		    ((std::uint32_t)(data[0])<<24) |
		    ((std::uint32_t)(data[1])<<16) |
		    ((std::uint32_t)(data[2])<<8) |
		    ((std::uint32_t)(data[3])<<0);
		std::memcpy(val, &v2, sizeof v2);
	}
}

static
void
parse(double* val, unsigned char const* const data, bool const littleEndian)
{
	assert(sizeof(double) == sizeof(std::uint64_t));
	if (littleEndian) {
		auto const v2 =
		    ((std::uint64_t)(data[0])<<0) |
		    ((std::uint64_t)(data[1])<<8) |
		    ((std::uint64_t)(data[2])<<16) |
		    ((std::uint64_t)(data[3])<<24) |
		    ((std::uint64_t)(data[4])<<32) |
		    ((std::uint64_t)(data[5])<<40) |
		    ((std::uint64_t)(data[6])<<48) |
		    ((std::uint64_t)(data[7])<<56);
		std::memcpy(val, &v2, sizeof v2);
	} else {
		auto const v2 =
		    ((std::uint64_t)(data[0])<<56) |
		    ((std::uint64_t)(data[1])<<48) |
		    ((std::uint64_t)(data[2])<<40) |
		    ((std::uint64_t)(data[3])<<32) |
		    ((std::uint64_t)(data[4])<<24) |
		    ((std::uint64_t)(data[5])<<16) |
		    ((std::uint64_t)(data[6])<<8) |
		    ((std::uint64_t)(data[7])<<0);
		std::memcpy(val, &v2, sizeof v2);
	}
}

static
int
parseHeader(
    Header* header,
    unsigned char const* const data,
    Epoch const epoch,
    Endian endian)
{
	std::int32_t nmag;
	if (endian == Endian::AUTO) {
		parse(&nmag, data + 0 + 4 + 4 + 4 + 4 + 4, true);
		if (std::abs(nmag) > 10) {
			parse(&nmag, data + 0 + 4 + 4 + 4 + 4 + 4, false);
			if (std::abs(nmag) > 10) {
				std::fprintf(stderr, "sidus: invalid header\n");
				return -1;
			} else {
				endian = Endian::BIG;
			}
		} else {
			endian = Endian::LITTLE;
		}
	}
	else if (endian == Endian::LITTLE) {
		parse(&nmag, data + 0 + 4 + 4 + 4 + 4 + 4, true);
		if (std::abs(nmag) > 10) {
				std::fprintf(stderr, "sidus: invalid header, maybe try big-endian?\n");
				return -1;
		}
	}
	else if (endian == Endian::BIG) {
		parse(&nmag, data + 0 + 4 + 4 + 4 + 4 + 4, false);
		if (std::abs(nmag) > 10) {
				std::fprintf(stderr, "sidus: invalid header, maybe try little-endian?\n");
				return -1;
		}
	}

	auto const littleEndian = endian == Endian::LITTLE;

	std::int32_t starn;
	parse(&starn, data + 0 + 4 + 4, littleEndian);
	std::int32_t stnum;
	parse(&stnum, data + 0 + 4 + 4 + 4, littleEndian);
	std::int32_t mprop;
	parse(&mprop, data + 0 + 4 + 4 + 4 + 4, littleEndian);
	std::int32_t nbent;
	parse(&nbent, data + 0 + 4 + 4 + 4 + 4 + 4 + 4, littleEndian);

	auto const isJ2000 = starn < 0 || nmag < 0;
	if (epoch == Epoch::J2000 && !isJ2000) {
		std::fprintf(stderr, "sidus: expected J2000 epoch but found B1950 epoch\n");
		return -1;
	}
	else if (epoch == Epoch::B1950 && isJ2000) {
		std::fprintf(stderr, "sidus: expected B1950 epoch but found J2000 epoch\n");
		return -1;
	}

	header->numStars = std::abs(starn);
	header->starId = stnum < 0 ? Header::NO_STAR_ID : (Header::StarId)stnum;
	header->starNameLength = stnum < 0 ? -stnum : 0;
	header->properMotion = (Header::ProperMotion)mprop;
	header->numMagnitudes = std::abs(nmag);
	header->numBytesPerStar = nbent;
	header->epoch = isJ2000 ? Epoch::J2000 : Epoch::B1950;
	header->littleEndian = littleEndian;

	return 0;
}

static
int
parseStar(
    Star* star,
    Header const& header,
    unsigned char const* const data)
{
	auto const littleEndian = header.littleEndian;
	auto cursor = 0;

	double xno = 0.0;
	if (header.starId == Header::CATALOG_STAR_ID ||
	    header.starId == Header::GSC_STAR_ID ||
	    header.starId == Header::TYCHO_STAR_ID) {
		float xno2;
		parse(&xno2, data + cursor, littleEndian);
		xno = xno2;
		cursor += 4;
	}
        else if (header.starId == Header::INTEGER_STAR_ID) {
		std::int32_t xno2;
		parse(&xno2, data + cursor, littleEndian);
		xno = xno2;
		cursor += 4;
	}

	double ra;
	parse(&ra, data + cursor, littleEndian);
	cursor += 8;
	double decl;
	parse(&decl, data + cursor, littleEndian);
	cursor += 8;
	char isp[2];
	isp[0] = *(data + cursor++);
	isp[1] = *(data + cursor++);
	std::int16_t mag;
	for (auto i = 0; i < header.numMagnitudes; ++i, cursor += 2) {
		if (i == header.apparentMagnitude) {
			parse(&mag, data + cursor, littleEndian);
		}
	}
	float xrpm = 0.0f;
	float xdpm = 0.0f;
	double svel = 0.0;
	if (header.properMotion == Header::PROPER_MOTION) {
		parse(&xrpm, data + cursor, littleEndian);
		cursor += 4;
		parse(&xdpm, data + cursor, littleEndian);
		cursor += 4;
	}
	else if (header.properMotion == Header::RADIAL_VELOCITY) {
		parse(&svel, data + cursor, littleEndian);
		cursor += 8;
	}
	char starname[header.starNameLength + 1];
	if (header.starNameLength > 0) {
		std::strncpy(starname, (char const*)(data + cursor), header.starNameLength);
	}
	starname[header.starNameLength] = '\0';

	star->name = starname;
	star->rightAscension = ra;
	star->declination = decl;
	star->starId = xno;
	star->magnitude = (float)(mag)/100.0f;
	star->properMotion.rightAscension = xrpm;
	star->properMotion.declination = xdpm;
	star->radialVelocity = svel;
	star->spectralType[0] = isp[0];
	star->spectralType[1] = isp[1];
	star->spectralType[2] = '\0';

	return 0;
}

static
std::string
sanitizeForC(char const * cs)
{
	std::string s;
	for (auto i = 0u; cs[i]; ++i) {
		auto const c = std::tolower(cs[i]);
		if (i == 0) {
			if (isalpha(c)) {
				s.push_back(c);
			} else {
				s.push_back('x');
			}
		} else {
			if (isalnum(c)) {
				s.push_back(c);
			} else {
				s.push_back('_');
			}
		}
	}
	return s;
}

static
void
printCHeader(char const* const inputfile,
	     unsigned const numStars,
	     Epoch const epoch,
	     bool const usefloat,
	     bool const usename,
	     bool const usetype)
{
	auto const var = sanitizeForC(inputfile);

	std::fprintf(stdout,
		     "/*\n"
		     " * Auto-generated from catalog %s by the sidus program\n"
		     " *\n"
		     " * Do this:\n"
		     " *   #define SIDUS_IMPLEMENTATION\n"
		     " * before you include this file in *one* C or C++ file to create the implementation\n"
		     " *\n"
		     " */\n\n",
		     inputfile);
	std::fprintf(stdout,
		     "#ifndef %s_h\n"
		     "#define %s_h\n\n",
		     var.c_str(), var.c_str());
	std::fputs("#ifdef __cplusplus\n"
		   "extern \"C\" {\n"
		   "#endif\n\n",
		   stdout);

	std::fputs("struct Star {\n", stdout);
	auto const epochstr = epoch == Epoch::J2000 ? "J2000" : "B1950";
	if (usefloat) {
		std::fprintf(stdout,
			     "	float rightAscension;	/* radians, %s */\n"
			     "	float declination;	/* radians, %s */\n"
			     "	float magnitude;\n",
			     epochstr, epochstr);
	} else {
		std::fprintf(stdout,
			     "	double rightAscension;	/* radians, %s */\n"
			     "	double declination;	/* radians, %s */\n"
			     "	double magnitude;\n",
			     epochstr, epochstr);
	}
	if (usename) {
		std::fputs("	const char *name;\n", stdout);
	}
	if (usetype) {
		std::fputs("	const char *type;\n", stdout);
	}
	std::fprintf(stdout, 
		     "};\n\n"
		     "enum { %s_num_stars = %u };\n\n"
		     "#ifndef SIDUS_IMPLEMENTATION\n"
		     "extern const struct Star * %s_stars;\n"
		     "#else\n"
		     "const struct Star %s_stars[%u] = {",
		     var.c_str(), numStars,
		     var.c_str(),
		     var.c_str(), numStars);
}

static
void
printCFooter()
{
	std::fputs("\n};\n\n"
		   "#endif\n\n"
		   "#ifdef __cplusplus\n"
		   "}\n"
		   "#endif\n\n"
		   "#endif\n",
		   stdout);
}

static
void
print(
    Star const& star,
    Header const & header,
    int const idx,
    bool const cformat,
    bool const usefloat,
    bool const usename,
    bool const usetype)
{
	if (cformat) {
		if (idx != 0) {
			std::fputs(", ", stdout);
		}
		if (usefloat) {
			std::fprintf(stdout,
				     "\n	{ % .9f, % .9f, % .9f",
				     star.rightAscension,
				     star.declination,
				     star.magnitude);
		} else {
			std::fprintf(stdout,
				     "\n	{ % .17lf, % .17lf, % .17lf",
				     star.rightAscension,
				     star.declination,
				     star.magnitude);
		}
		if (usename) {
			std::fprintf(stdout,
				     ", \"%s\"",
				     star.name.c_str());
		}
		if (usetype) {
			std::fprintf(stdout,
				     ", \"%s\"",
				     star.spectralType);
		}
		std::fputs(" }", stdout);
	} else {
		if (usename) {
			std::fprintf(stdout, "%s,", star.name.c_str());
		}
		if (usefloat) {
			std::fprintf(stdout,
				     "%.9f,%.9f,%.9f",
				     star.rightAscension,
				     star.declination,
				     star.magnitude);
		} else {
			std::fprintf(stdout,
				     "%.17lf,%.17lf,%.17lf",
				     star.rightAscension,
				     star.declination,
				     star.magnitude);
		}
		if (usetype) {
			std::fprintf(stdout,
				     ",%c%c",
				     star.spectralType[0],
				     star.spectralType[1]);
		}
		std::fputc('\n', stdout);
	}
}

}	// !namespace

int
main(int argc, char** argv)
{
	if (argc < 1) {
		usage(stderr);
		return -1;
	}

	auto apparentMagnitude = 0;
	auto filterMagnitude = DBL_MAX;
	Epoch epoch = Epoch::AUTO;
	auto cformat = false;
	Endian endian = Endian::AUTO;
	auto usefloat = false;
	auto onlymeta = false;
	enum class Sort { NO, MAG, RA } sort = Sort::NO;
	auto usename = false;
	auto usetype = false;
	char const* inputfile = nullptr;

	for (auto i = 1; i < argc; ++i) {
		auto const arg = std::string(argv[i]);
		if (arg[0] == '-') {
			if (arg.size() < 2) {
				usage(stderr);
				return -1;
			}
			switch (arg[1]) {
			case 'a':
				if (arg.size() < 3) {
					std::fprintf(stderr, "Invalid option '%s'\n", arg.c_str());
					usage(stderr);
					return -1;
				}
				apparentMagnitude = std::stod(arg.substr(2));
				continue;
			case 'f':
				if (arg.size() < 3) {
					std::fprintf(stderr, "Invalid option '%s'\n", arg.c_str());
					usage(stderr);
					return -1;
				}
				filterMagnitude = std::stod(arg.substr(2));
				continue;
			case 'B':
				if (arg.size() < 3 || arg.substr(2) != "1950") {
					std::fprintf(stderr, "Invalid option '%s'\n", arg.c_str());
					usage(stderr);
					return -1;
				}
				epoch = Epoch::B1950;
				continue;
			case 'J':
				if (arg.size() < 3 || arg.substr(2) != "2000") {
					std::fprintf(stderr, "Invalid option '%s'\n", arg.c_str());
					usage(stderr);
					return -1;
				}
				epoch = Epoch::J2000;
				continue;
			case 'c':
				cformat = true;
				break;
			case 'l':
				if (arg.size() < 3 || arg[2] != 'e') {
					std::fprintf(stderr, "Invalid option '%s'\n", arg.c_str());
					usage(stderr);
					return -1;
				}
				endian = Endian::LITTLE;
				continue;
			case 'b':
				if (arg.size() < 3 || arg[2] != 'e') {
					std::fprintf(stderr, "Invalid option '%s'\n", arg.c_str());
					usage(stderr);
					return -1;
				}
				endian = Endian::BIG;
				continue;
			case 's':
				usefloat = true;
				break;
			case 'i':
				onlymeta = true;
				break;
			case 'm':
				sort = Sort::MAG;
				break;
			case 'r':
				sort = Sort::RA;
				break;
			case 'n':
				usename = true;
				break;
			case 'p':
				usetype = true;
				break;
			case 'h':
				usage(stdout);
				return 0;
			case 'v':
				version();
				return 0;
			case '-':
				{
					auto const larg = arg.substr(2);
					if (larg == "help") {
						usage(stdout);
						return 0;
					}
					else if (larg == "version") {
						version();
						return 0;

					} else {
						usage(stderr);
						return -1;
					}
				}
				continue;
			default:
				usage(stderr);
				return -1;
			}
			if (arg.size() != 2) {
				std::fprintf(stderr, "Invalid option '%s'\n", arg.c_str());
				usage(stderr);
				return -1;
			}
		} else {
			inputfile = argv[i];
		}
	}

	if (!inputfile) {
		std::fprintf(stderr, "sidus: no input file\n");
		usage(stderr);
		return -1;
	}

	size_t filesize = 0;
	if (getFileSize(&filesize, inputfile) != 0) {
		std::fprintf(stderr, "sidus: %s: failed to open file\n", inputfile);
		return -1;
	}

	if (filesize == 0) {
		std::fprintf(stderr, "sidus: %s: empty\n", inputfile);
		return -1;
	}

	if (filesize < 28) {
		std::fprintf(stderr, "sidus: %s: no header\n", inputfile);
		return -1;
	}


	auto data = std::unique_ptr<unsigned char[]>(new unsigned char[filesize]);
	if (readFile(data.get(), inputfile, filesize) != 0) {
		std::fprintf(stderr, "sidus: %s: failed to read file\n", inputfile);
		return -1;
	}

	Header header;
	if (parseHeader(&header, data.get(), epoch, endian) != 0) {
		return -1;
	}

	auto const starDataSize = header.numStars*header.numBytesPerStar;

	if (filesize < (28 + starDataSize)) {
		std::fprintf(stderr, "sidus: header.numStars: %d, bytesPerStar: %d, %lu < %d, file too short\n",
			     header.numStars, header.numBytesPerStar, filesize, 28 + starDataSize);
		return -1;
	}
	if (header.numMagnitudes < 1) {
		std::fprintf(stderr, "sidus: expected at least one magnitude per star, found: %d\n",
			     header.numMagnitudes);
		return -1;
	}
	header.apparentMagnitude = std::min(header.numMagnitudes - 1, header.numMagnitudes);
	if (header.starNameLength == 0) {
		usename = false;
	}

	if (onlymeta) {
		fprintf(stdout,
			"Catalog information:\n"
			" Number of stars: %d\n"
			" Id: %s\n"
			" Names: %s\n"
			" Proper motion: %s\n"
			" Number of magnitudes: %d\n"
			" Epoch: %s\n"
			" Bytes per star: %d\n",
			header.numStars,
			header.starId == Header::NO_STAR_ID ?
				"No" : header.starId == Header::CATALOG_STAR_ID ?
				"Catalog star id" : header.starId == Header::GSC_STAR_ID ?
				"GSC star id" : header.starId == Header::TYCHO_STAR_ID ?
				"Tycho star id" : header.starId == Header::INTEGER_STAR_ID ?
				"Integer star id" : "UNKNOWN",
			header.starNameLength > 0 ? "Yes" : "No",
			header.properMotion == Header::NO_PROPER_MOTION ?
				"No" : header.properMotion == Header::PROPER_MOTION ?
				"Yes" : header.properMotion == Header::RADIAL_VELOCITY ?
				"Radial velocity" : "UNKNOWN",
			header.numMagnitudes,
			header.epoch == Epoch::J2000 ? "J2000" : "B1950",
			header.numBytesPerStar
			);
		return 0;
	}

	auto cursor = 28;
	std::multimap<double, Star> map;
	for (auto i = 0; i < header.numStars; ++i, cursor += header.numBytesPerStar) {
		Star star;
		if (parseStar(&star, header, data.get() + cursor) != 0) {
			continue;
		}
		if (star.magnitude > filterMagnitude) {
			continue;
		}
		// Filter out "invalid" entries
		if (star.magnitude == 0.0 &&
		    star.rightAscension == 0.0 &&
		    star.declination == 0.0) {
			continue;
		}

		switch (sort) {
		case Sort::RA:
			map.insert(std::make_pair(star.rightAscension, star));
			break;
		case Sort::MAG:
			map.insert(std::make_pair(star.magnitude, star));
			break;
		default:
			map.insert(std::make_pair(i, star));
			break;
		}
	}

	if (cformat) {
		printCHeader(inputfile, map.size(), header.epoch, usefloat, usename, usetype);
	}

	auto idx = 0;
	for (auto const & pair : map) {
		print(pair.second, header, idx++, cformat, usefloat, usename, usetype);
	}

	if (cformat) {
		printCFooter();
	}
}

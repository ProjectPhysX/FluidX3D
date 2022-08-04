#pragma once

/*
LodePNG version 20200306

Copyright (c) 2005-2020 Lode Vandevenne

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

	1. The origin of this software must not be misrepresented; you must not
	claim that you wrote the original software. If you use this software
	in a product, an acknowledgment in the product documentation would be
	appreciated but is not required.

	2. Altered source versions must be plainly marked as such, and must not be
	misrepresented as being the original software.

	3. This notice may not be removed or altered from any source
	distribution.
*/

#ifndef LODEPNG_H
#define LODEPNG_H

#include <string.h> /*for size_t*/

extern const char* LODEPNG_VERSION_STRING;

/*
The following #defines are used to create code sections. They can be disabled
to disable code sections, which can give faster compile time and smaller binary.
The "NO_COMPILE" defines are designed to be used to pass as defines to the
compiler command to disable them without modifying this header, e.g.
-DLODEPNG_NO_COMPILE_ZLIB for gcc.
In addition to those below, you can also define LODEPNG_NO_COMPILE_CRC to
allow implementing a custom lodepng_crc32.
*/
/*deflate & zlib. If disabled, you must specify alternative zlib functions in
the custom_zlib field of the compress and decompress settings*/
#ifndef LODEPNG_NO_COMPILE_ZLIB
#define LODEPNG_COMPILE_ZLIB
#endif

/*png encoder and png decoder*/
#ifndef LODEPNG_NO_COMPILE_PNG
#define LODEPNG_COMPILE_PNG
#endif

/*deflate&zlib decoder and png decoder*/
#ifndef LODEPNG_NO_COMPILE_DECODER
#define LODEPNG_COMPILE_DECODER
#endif

/*deflate&zlib encoder and png encoder*/
#ifndef LODEPNG_NO_COMPILE_ENCODER
#define LODEPNG_COMPILE_ENCODER
#endif

/*the optional built in harddisk file loading and saving functions*/
#ifndef LODEPNG_NO_COMPILE_DISK
#define LODEPNG_COMPILE_DISK
#endif

/*support for chunks other than IHDR, IDAT, PLTE, tRNS, IEND: ancillary and unknown chunks*/
#ifndef LODEPNG_NO_COMPILE_ANCILLARY_CHUNKS
#define LODEPNG_COMPILE_ANCILLARY_CHUNKS
#endif

/*ability to convert error numerical codes to English text string*/
#ifndef LODEPNG_NO_COMPILE_ERROR_TEXT
#define LODEPNG_COMPILE_ERROR_TEXT
#endif

/*Compile the default allocators (C's free, malloc and realloc). If you disable this,
you can define the functions lodepng_free, lodepng_malloc and lodepng_realloc in your
source files with custom allocators.*/
#ifndef LODEPNG_NO_COMPILE_ALLOCATORS
#define LODEPNG_COMPILE_ALLOCATORS
#endif

/*compile the C++ version (you can disable the C++ wrapper here even when compiling for C++)*/
#ifdef __cplusplus
#ifndef LODEPNG_NO_COMPILE_CPP
#define LODEPNG_COMPILE_CPP
#endif
#endif

#ifdef LODEPNG_COMPILE_CPP
#include <vector>
#include <string>
#endif /*LODEPNG_COMPILE_CPP*/

#ifdef LODEPNG_COMPILE_PNG
/*The PNG color types (also used for raw image).*/
typedef enum LodePNGColorType {
	LCT_GREY = 0, /*grayscale: 1,2,4,8,16 bit*/
	LCT_RGB = 2, /*RGB: 8,16 bit*/
	LCT_PALETTE = 3, /*palette: 1,2,4,8 bit*/
	LCT_GREY_ALPHA = 4, /*grayscale with alpha: 8,16 bit*/
	LCT_RGBA = 6, /*RGB with alpha: 8,16 bit*/
	/*LCT_MAX_OCTET_VALUE lets the compiler allow this enum to represent any invalid
	byte value from 0 to 255 that could be present in an invalid PNG file header. Do
	not use, compare with or set the name LCT_MAX_OCTET_VALUE, instead either use
	the valid color type names above, or numeric values like 1 or 7 when checking for
	particular disallowed color type byte values, or cast to integer to print it.*/
	LCT_MAX_OCTET_VALUE = 255
} LodePNGColorType;

#ifdef LODEPNG_COMPILE_DECODER
/*
Converts PNG data in memory to raw pixel data.
out: Output parameter. Pointer to buffer that will contain the raw pixel data.
	 After decoding, its size is w * h * (bytes per pixel) bytes larger than
	 initially. Bytes per pixel depends on colortype and bitdepth.
	 Must be freed after usage with free(*out).
	 Note: for 16-bit per channel colors, uses big endian format like PNG does.
w: Output parameter. Pointer to width of pixel data.
h: Output parameter. Pointer to height of pixel data.
in: Memory buffer with the PNG file.
insize: size of the in buffer.
colortype: the desired color type for the raw output image. See explanation on PNG color types.
bitdepth: the desired bit depth for the raw output image. See explanation on PNG color types.
Return value: LodePNG error code (0 means no error).
*/
unsigned lodepng_decode_memory(unsigned char** out, unsigned* w, unsigned* h,
	const unsigned char* in, size_t insize,
	LodePNGColorType colortype, unsigned bitdepth);

/*Same as lodepng_decode_memory, but always decodes to 32-bit RGBA raw image*/
unsigned lodepng_decode32(unsigned char** out, unsigned* w, unsigned* h,
	const unsigned char* in, size_t insize);

/*Same as lodepng_decode_memory, but always decodes to 24-bit RGB raw image*/
unsigned lodepng_decode24(unsigned char** out, unsigned* w, unsigned* h,
	const unsigned char* in, size_t insize);

#ifdef LODEPNG_COMPILE_DISK
/*
Load PNG from disk, from file with given name.
Same as the other decode functions, but instead takes a filename as input.
*/
unsigned lodepng_decode_file(unsigned char** out, unsigned* w, unsigned* h,
	const char* filename,
	LodePNGColorType colortype, unsigned bitdepth);

/*Same as lodepng_decode_file, but always decodes to 32-bit RGBA raw image.*/
unsigned lodepng_decode32_file(unsigned char** out, unsigned* w, unsigned* h,
	const char* filename);

/*Same as lodepng_decode_file, but always decodes to 24-bit RGB raw image.*/
unsigned lodepng_decode24_file(unsigned char** out, unsigned* w, unsigned* h,
	const char* filename);
#endif /*LODEPNG_COMPILE_DISK*/
#endif /*LODEPNG_COMPILE_DECODER*/


#ifdef LODEPNG_COMPILE_ENCODER
/*
Converts raw pixel data into a PNG image in memory. The colortype and bitdepth
  of the output PNG image cannot be chosen, they are automatically determined
  by the colortype, bitdepth and content of the input pixel data.
  Note: for 16-bit per channel colors, needs big endian format like PNG does.
out: Output parameter. Pointer to buffer that will contain the PNG image data.
	 Must be freed after usage with free(*out).
outsize: Output parameter. Pointer to the size in bytes of the out buffer.
image: The raw pixel data to encode. The size of this buffer should be
	   w * h * (bytes per pixel), bytes per pixel depends on colortype and bitdepth.
w: width of the raw pixel data in pixels.
h: height of the raw pixel data in pixels.
colortype: the color type of the raw input image. See explanation on PNG color types.
bitdepth: the bit depth of the raw input image. See explanation on PNG color types.
Return value: LodePNG error code (0 means no error).
*/
unsigned lodepng_encode_memory(unsigned char** out, size_t* outsize,
	const unsigned char* image, unsigned w, unsigned h,
	LodePNGColorType colortype, unsigned bitdepth);

/*Same as lodepng_encode_memory, but always encodes from 32-bit RGBA raw image.*/
unsigned lodepng_encode32(unsigned char** out, size_t* outsize,
	const unsigned char* image, unsigned w, unsigned h);

/*Same as lodepng_encode_memory, but always encodes from 24-bit RGB raw image.*/
unsigned lodepng_encode24(unsigned char** out, size_t* outsize,
	const unsigned char* image, unsigned w, unsigned h);

#ifdef LODEPNG_COMPILE_DISK
/*
Converts raw pixel data into a PNG file on disk.
Same as the other encode functions, but instead takes a filename as output.
NOTE: This overwrites existing files without warning!
*/
unsigned lodepng_encode_file(const char* filename,
	const unsigned char* image, unsigned w, unsigned h,
	LodePNGColorType colortype, unsigned bitdepth);

/*Same as lodepng_encode_file, but always encodes from 32-bit RGBA raw image.*/
unsigned lodepng_encode32_file(const char* filename,
	const unsigned char* image, unsigned w, unsigned h);

/*Same as lodepng_encode_file, but always encodes from 24-bit RGB raw image.*/
unsigned lodepng_encode24_file(const char* filename,
	const unsigned char* image, unsigned w, unsigned h);
#endif /*LODEPNG_COMPILE_DISK*/
#endif /*LODEPNG_COMPILE_ENCODER*/


#ifdef LODEPNG_COMPILE_CPP
namespace lodepng {
#ifdef LODEPNG_COMPILE_DECODER
	/*Same as lodepng_decode_memory, but decodes to an std::vector. The colortype
	is the format to output the pixels to. Default is RGBA 8-bit per channel.*/
	unsigned decode(std::vector<unsigned char>& out, unsigned& w, unsigned& h,
		const unsigned char* in, size_t insize,
		LodePNGColorType colortype = LCT_RGBA, unsigned bitdepth = 8);
	unsigned decode(std::vector<unsigned char>& out, unsigned& w, unsigned& h,
		const std::vector<unsigned char>& in,
		LodePNGColorType colortype = LCT_RGBA, unsigned bitdepth = 8);
#ifdef LODEPNG_COMPILE_DISK
	/*
	Converts PNG file from disk to raw pixel data in memory.
	Same as the other decode functions, but instead takes a filename as input.
	*/
	unsigned decode(std::vector<unsigned char>& out, unsigned& w, unsigned& h,
		const std::string& filename,
		LodePNGColorType colortype = LCT_RGBA, unsigned bitdepth = 8);
#endif /* LODEPNG_COMPILE_DISK */
#endif /* LODEPNG_COMPILE_DECODER */

#ifdef LODEPNG_COMPILE_ENCODER
	/*Same as lodepng_encode_memory, but encodes to an std::vector. colortype
	is that of the raw input data. The output PNG color type will be auto chosen.*/
	unsigned encode(std::vector<unsigned char>& out,
		const unsigned char* in, unsigned w, unsigned h,
		LodePNGColorType colortype = LCT_RGBA, unsigned bitdepth = 8);
	unsigned encode(std::vector<unsigned char>& out,
		const std::vector<unsigned char>& in, unsigned w, unsigned h,
		LodePNGColorType colortype = LCT_RGBA, unsigned bitdepth = 8);
#ifdef LODEPNG_COMPILE_DISK
	/*
	Converts 32-bit RGBA raw pixel data into a PNG file on disk.
	Same as the other encode functions, but instead takes a filename as output.
	NOTE: This overwrites existing files without warning!
	*/
	unsigned encode(const std::string& filename,
		const unsigned char* in, unsigned w, unsigned h,
		LodePNGColorType colortype = LCT_RGBA, unsigned bitdepth = 8);
	unsigned encode(const std::string& filename,
		const std::vector<unsigned char>& in, unsigned w, unsigned h,
		LodePNGColorType colortype = LCT_RGBA, unsigned bitdepth = 8);
#endif /* LODEPNG_COMPILE_DISK */
#endif /* LODEPNG_COMPILE_ENCODER */
} /* namespace lodepng */
#endif /*LODEPNG_COMPILE_CPP*/
#endif /*LODEPNG_COMPILE_PNG*/

#ifdef LODEPNG_COMPILE_ERROR_TEXT
/*Returns an English description of the numerical error code.*/
const char* lodepng_error_text(unsigned code);
#endif /*LODEPNG_COMPILE_ERROR_TEXT*/

#ifdef LODEPNG_COMPILE_DECODER
/*Settings for zlib decompression*/
typedef struct LodePNGDecompressSettings LodePNGDecompressSettings;
struct LodePNGDecompressSettings {
	/* Check LodePNGDecoderSettings for more ignorable errors such as ignore_crc */
	unsigned ignore_adler32; /*if 1, continue and don't give an error message if the Adler32 checksum is corrupted*/
	unsigned ignore_nlen; /*ignore complement of len checksum in uncompressed blocks*/

	/*use custom zlib decoder instead of built in one (default: null)*/
	unsigned(*custom_zlib)(unsigned char**, size_t*,
		const unsigned char*, size_t,
		const LodePNGDecompressSettings*);
	/*use custom deflate decoder instead of built in one (default: null)
	if custom_zlib is not null, custom_inflate is ignored (the zlib format uses deflate)*/
	unsigned(*custom_inflate)(unsigned char**, size_t*,
		const unsigned char*, size_t,
		const LodePNGDecompressSettings*);

	const void* custom_context; /*optional custom settings for custom functions*/
};

extern const LodePNGDecompressSettings lodepng_default_decompress_settings;
void lodepng_decompress_settings_init(LodePNGDecompressSettings* settings);
#endif /*LODEPNG_COMPILE_DECODER*/

#ifdef LODEPNG_COMPILE_ENCODER
/*
Settings for zlib compression. Tweaking these settings tweaks the balance
between speed and compression ratio.
*/
typedef struct LodePNGCompressSettings LodePNGCompressSettings;
struct LodePNGCompressSettings /*deflate = compress*/ {
	/*LZ77 related settings*/
	unsigned btype; /*the block type for LZ (0, 1, 2 or 3, see zlib standard). Should be 2 for proper compression.*/
	unsigned use_lz77; /*whether or not to use LZ77. Should be 1 for proper compression.*/
	unsigned windowsize; /*must be a power of two <= 32768. higher compresses more but is slower. Default value: 2048.*/
	unsigned minmatch; /*minimum lz77 length. 3 is normally best, 6 can be better for some PNGs. Default: 0*/
	unsigned nicematch; /*stop searching if >= this length found. Set to 258 for best compression. Default: 128*/
	unsigned lazymatching; /*use lazy matching: better compression but a bit slower. Default: true*/

	/*use custom zlib encoder instead of built in one (default: null)*/
	unsigned(*custom_zlib)(unsigned char**, size_t*,
		const unsigned char*, size_t,
		const LodePNGCompressSettings*);
	/*use custom deflate encoder instead of built in one (default: null)
	if custom_zlib is used, custom_deflate is ignored since only the built in
	zlib function will call custom_deflate*/
	unsigned(*custom_deflate)(unsigned char**, size_t*,
		const unsigned char*, size_t,
		const LodePNGCompressSettings*);

	const void* custom_context; /*optional custom settings for custom functions*/
};

extern const LodePNGCompressSettings lodepng_default_compress_settings;
void lodepng_compress_settings_init(LodePNGCompressSettings* settings);
#endif /*LODEPNG_COMPILE_ENCODER*/

#ifdef LODEPNG_COMPILE_PNG
/*
Color mode of an image. Contains all information required to decode the pixel
bits to RGBA colors. This information is the same as used in the PNG file
format, and is used both for PNG and raw image data in LodePNG.
*/
typedef struct LodePNGColorMode {
	/*header (IHDR)*/
	LodePNGColorType colortype; /*color type, see PNG standard or documentation further in this header file*/
	unsigned bitdepth;  /*bits per sample, see PNG standard or documentation further in this header file*/

	/*
	palette (PLTE and tRNS)

	Dynamically allocated with the colors of the palette, including alpha.
	This field may not be allocated directly, use lodepng_color_mode_init first,
	then lodepng_palette_add per color to correctly initialize it (to ensure size
	of exactly 1024 bytes).

	The alpha channels must be set as well, set them to 255 for opaque images.

	When decoding, by default you can ignore this palette, since LodePNG already
	fills the palette colors in the pixels of the raw RGBA output.

	The palette is only supported for color type 3.
	*/
	unsigned char* palette; /*palette in RGBARGBA... order. Must be either 0, or when allocated must have 1024 bytes*/
	size_t palettesize; /*palette size in number of colors (amount of used bytes is 4 * palettesize)*/

	/*
	transparent color key (tRNS)

	This color uses the same bit depth as the bitdepth value in this struct, which can be 1-bit to 16-bit.
	For grayscale PNGs, r, g and b will all 3 be set to the same.

	When decoding, by default you can ignore this information, since LodePNG sets
	pixels with this key to transparent already in the raw RGBA output.

	The color key is only supported for color types 0 and 2.
	*/
	unsigned key_defined; /*is a transparent color key given? 0 = false, 1 = true*/
	unsigned key_r;       /*red/grayscale component of color key*/
	unsigned key_g;       /*green component of color key*/
	unsigned key_b;       /*blue component of color key*/
} LodePNGColorMode;

/*init, cleanup and copy functions to use with this struct*/
void lodepng_color_mode_init(LodePNGColorMode* info);
void lodepng_color_mode_cleanup(LodePNGColorMode* info);
/*return value is error code (0 means no error)*/
unsigned lodepng_color_mode_copy(LodePNGColorMode* dest, const LodePNGColorMode* source);
/* Makes a temporary LodePNGColorMode that does not need cleanup (no palette) */
LodePNGColorMode lodepng_color_mode_make(LodePNGColorType colortype, unsigned bitdepth);

void lodepng_palette_clear(LodePNGColorMode* info);
/*add 1 color to the palette*/
unsigned lodepng_palette_add(LodePNGColorMode* info,
	unsigned char r, unsigned char g, unsigned char b, unsigned char a);

/*get the total amount of bits per pixel, based on colortype and bitdepth in the struct*/
unsigned lodepng_get_bpp(const LodePNGColorMode* info);
/*get the amount of color channels used, based on colortype in the struct.
If a palette is used, it counts as 1 channel.*/
unsigned lodepng_get_channels(const LodePNGColorMode* info);
/*is it a grayscale type? (only colortype 0 or 4)*/
unsigned lodepng_is_greyscale_type(const LodePNGColorMode* info);
/*has it got an alpha channel? (only colortype 2 or 6)*/
unsigned lodepng_is_alpha_type(const LodePNGColorMode* info);
/*has it got a palette? (only colortype 3)*/
unsigned lodepng_is_palette_type(const LodePNGColorMode* info);
/*only returns true if there is a palette and there is a value in the palette with alpha < 255.
Loops through the palette to check this.*/
unsigned lodepng_has_palette_alpha(const LodePNGColorMode* info);
/*
Check if the given color info indicates the possibility of having non-opaque pixels in the PNG image.
Returns true if the image can have translucent or invisible pixels (it still be opaque if it doesn't use such pixels).
Returns false if the image can only have opaque pixels.
In detail, it returns true only if it's a color type with alpha, or has a palette with non-opaque values,
or if "key_defined" is true.
*/
unsigned lodepng_can_have_alpha(const LodePNGColorMode* info);
/*Returns the byte size of a raw image buffer with given width, height and color mode*/
size_t lodepng_get_raw_size(unsigned w, unsigned h, const LodePNGColorMode* color);

#ifdef LODEPNG_COMPILE_ANCILLARY_CHUNKS
/*The information of a Time chunk in PNG.*/
typedef struct LodePNGTime {
	unsigned year;    /*2 bytes used (0-65535)*/
	unsigned month;   /*1-12*/
	unsigned day;     /*1-31*/
	unsigned hour;    /*0-23*/
	unsigned minute;  /*0-59*/
	unsigned second;  /*0-60 (to allow for leap seconds)*/
} LodePNGTime;
#endif /*LODEPNG_COMPILE_ANCILLARY_CHUNKS*/

/*Information about the PNG image, except pixels, width and height.*/
typedef struct LodePNGInfo {
	/*header (IHDR), palette (PLTE) and transparency (tRNS) chunks*/
	unsigned compression_method;/*compression method of the original file. Always 0.*/
	unsigned filter_method;     /*filter method of the original file*/
	unsigned interlace_method;  /*interlace method of the original file: 0=none, 1=Adam7*/
	LodePNGColorMode color;     /*color type and bits, palette and transparency of the PNG file*/

#ifdef LODEPNG_COMPILE_ANCILLARY_CHUNKS
  /*
  Suggested background color chunk (bKGD)

  This uses the same color mode and bit depth as the PNG (except no alpha channel),
  with values truncated to the bit depth in the unsigned integer.

  For grayscale and palette PNGs, the value is stored in background_r. The values
  in background_g and background_b are then unused.

  So when decoding, you may get these in a different color mode than the one you requested
  for the raw pixels.

  When encoding with auto_convert, you must use the color model defined in info_png.color for
  these values. The encoder normally ignores info_png.color when auto_convert is on, but will
  use it to interpret these values (and convert copies of them to its chosen color model).

  When encoding, avoid setting this to an expensive color, such as a non-gray value
  when the image is gray, or the compression will be worse since it will be forced to
  write the PNG with a more expensive color mode (when auto_convert is on).

  The decoder does not use this background color to edit the color of pixels. This is a
  completely optional metadata feature.
  */
	unsigned background_defined; /*is a suggested background color given?*/
	unsigned background_r;       /*red/gray/palette component of suggested background color*/
	unsigned background_g;       /*green component of suggested background color*/
	unsigned background_b;       /*blue component of suggested background color*/

	/*
	non-international text chunks (tEXt and zTXt)

	The char** arrays each contain num strings. The actual messages are in
	text_strings, while text_keys are keywords that give a short description what
	the actual text represents, e.g. Title, Author, Description, or anything else.

	All the string fields below including keys, names and language tags are null terminated.
	The PNG specification uses null characters for the keys, names and tags, and forbids null
	characters to appear in the main text which is why we can use null termination everywhere here.

	A keyword is minimum 1 character and maximum 79 characters long. It's
	discouraged to use a single line length longer than 79 characters for texts.

	Don't allocate these text buffers yourself. Use the init/cleanup functions
	correctly and use lodepng_add_text and lodepng_clear_text.
	*/
	size_t text_num; /*the amount of texts in these char** buffers (there may be more texts in itext)*/
	char** text_keys; /*the keyword of a text chunk (e.g. "Comment")*/
	char** text_strings; /*the actual text*/

	/*
	international text chunks (iTXt)
	Similar to the non-international text chunks, but with additional strings
	"langtags" and "transkeys".
	*/
	size_t itext_num; /*the amount of international texts in this PNG*/
	char** itext_keys; /*the English keyword of the text chunk (e.g. "Comment")*/
	char** itext_langtags; /*language tag for this text's language, ISO/IEC 646 string, e.g. ISO 639 language tag*/
	char** itext_transkeys; /*keyword translated to the international language - UTF-8 string*/
	char** itext_strings; /*the actual international text - UTF-8 string*/

	/*time chunk (tIME)*/
	unsigned time_defined; /*set to 1 to make the encoder generate a tIME chunk*/
	LodePNGTime time;

	/*phys chunk (pHYs)*/
	unsigned phys_defined; /*if 0, there is no pHYs chunk and the values below are undefined, if 1 else there is one*/
	unsigned phys_x; /*pixels per unit in x direction*/
	unsigned phys_y; /*pixels per unit in y direction*/
	unsigned phys_unit; /*may be 0 (unknown unit) or 1 (metre)*/

	/*
	Color profile related chunks: gAMA, cHRM, sRGB, iCPP

	LodePNG does not apply any color conversions on pixels in the encoder or decoder and does not interpret these color
	profile values. It merely passes on the information. If you wish to use color profiles and convert colors, please
	use these values with a color management library.

	See the PNG, ICC and sRGB specifications for more information about the meaning of these values.
	*/

	/* gAMA chunk: optional, overridden by sRGB or iCCP if those are present. */
	unsigned gama_defined; /* Whether a gAMA chunk is present (0 = not present, 1 = present). */
	unsigned gama_gamma;   /* Gamma exponent times 100000 */

	/* cHRM chunk: optional, overridden by sRGB or iCCP if those are present. */
	unsigned chrm_defined; /* Whether a cHRM chunk is present (0 = not present, 1 = present). */
	unsigned chrm_white_x; /* White Point x times 100000 */
	unsigned chrm_white_y; /* White Point y times 100000 */
	unsigned chrm_red_x;   /* Red x times 100000 */
	unsigned chrm_red_y;   /* Red y times 100000 */
	unsigned chrm_green_x; /* Green x times 100000 */
	unsigned chrm_green_y; /* Green y times 100000 */
	unsigned chrm_blue_x;  /* Blue x times 100000 */
	unsigned chrm_blue_y;  /* Blue y times 100000 */

	/*
	sRGB chunk: optional. May not appear at the same time as iCCP.
	If gAMA is also present gAMA must contain value 45455.
	If cHRM is also present cHRM must contain respectively 31270,32900,64000,33000,30000,60000,15000,6000.
	*/
	unsigned srgb_defined; /* Whether an sRGB chunk is present (0 = not present, 1 = present). */
	unsigned srgb_intent;  /* Rendering intent: 0=perceptual, 1=rel. colorimetric, 2=saturation, 3=abs. colorimetric */

	/*
	iCCP chunk: optional. May not appear at the same time as sRGB.

	LodePNG does not parse or use the ICC profile (except its color space header field for an edge case), a
	separate library to handle the ICC data (not included in LodePNG) format is needed to use it for color
	management and conversions.

	For encoding, if iCCP is present, gAMA and cHRM are recommended to be added as well with values that match the ICC
	profile as closely as possible, if you wish to do this you should provide the correct values for gAMA and cHRM and
	enable their '_defined' flags since LodePNG will not automatically compute them from the ICC profile.

	For encoding, the ICC profile is required by the PNG specification to be an "RGB" profile for non-gray
	PNG color types and a "GRAY" profile for gray PNG color types. If you disable auto_convert, you must ensure
	the ICC profile type matches your requested color type, else the encoder gives an error. If auto_convert is
	enabled (the default), and the ICC profile is not a good match for the pixel data, this will result in an encoder
	error if the pixel data has non-gray pixels for a GRAY profile, or a silent less-optimal compression of the pixel
	data if the pixels could be encoded as grayscale but the ICC profile is RGB.

	To avoid this do not set an ICC profile in the image unless there is a good reason for it, and when doing so
	make sure you compute it carefully to avoid the above problems.
	*/
	unsigned iccp_defined;      /* Whether an iCCP chunk is present (0 = not present, 1 = present). */
	char* iccp_name;            /* Null terminated string with profile name, 1-79 bytes */
	/*
	The ICC profile in iccp_profile_size bytes.
	Don't allocate this buffer yourself. Use the init/cleanup functions
	correctly and use lodepng_set_icc and lodepng_clear_icc.
	*/
	unsigned char* iccp_profile;
	unsigned iccp_profile_size; /* The size of iccp_profile in bytes */

	/* End of color profile related chunks */


	/*
	unknown chunks: chunks not known by LodePNG, passed on byte for byte.

	There are 3 buffers, one for each position in the PNG where unknown chunks can appear.
	Each buffer contains all unknown chunks for that position consecutively.
	The 3 positions are:
	0: between IHDR and PLTE, 1: between PLTE and IDAT, 2: between IDAT and IEND.

	For encoding, do not store critical chunks or known chunks that are enabled with a "_defined" flag
	above in here, since the encoder will blindly follow this and could then encode an invalid PNG file
	(such as one with two IHDR chunks or the disallowed combination of sRGB with iCCP). But do use
	this if you wish to store an ancillary chunk that is not supported by LodePNG (such as sPLT or hIST),
	or any non-standard PNG chunk.

	Do not allocate or traverse this data yourself. Use the chunk traversing functions declared
	later, such as lodepng_chunk_next and lodepng_chunk_append, to read/write this struct.
	*/
	unsigned char* unknown_chunks_data[3];
	size_t unknown_chunks_size[3]; /*size in bytes of the unknown chunks, given for protection*/
#endif /*LODEPNG_COMPILE_ANCILLARY_CHUNKS*/
} LodePNGInfo;

/*init, cleanup and copy functions to use with this struct*/
void lodepng_info_init(LodePNGInfo* info);
void lodepng_info_cleanup(LodePNGInfo* info);
/*return value is error code (0 means no error)*/
unsigned lodepng_info_copy(LodePNGInfo* dest, const LodePNGInfo* source);

#ifdef LODEPNG_COMPILE_ANCILLARY_CHUNKS
unsigned lodepng_add_text(LodePNGInfo* info, const char* key, const char* str); /*push back both texts at once*/
void lodepng_clear_text(LodePNGInfo* info); /*use this to clear the texts again after you filled them in*/

unsigned lodepng_add_itext(LodePNGInfo* info, const char* key, const char* langtag,
	const char* transkey, const char* str); /*push back the 4 texts of 1 chunk at once*/
void lodepng_clear_itext(LodePNGInfo* info); /*use this to clear the itexts again after you filled them in*/

/*replaces if exists*/
unsigned lodepng_set_icc(LodePNGInfo* info, const char* name, const unsigned char* profile, unsigned profile_size);
void lodepng_clear_icc(LodePNGInfo* info); /*use this to clear the texts again after you filled them in*/
#endif /*LODEPNG_COMPILE_ANCILLARY_CHUNKS*/

/*
Converts raw buffer from one color type to another color type, based on
LodePNGColorMode structs to describe the input and output color type.
See the reference manual at the end of this header file to see which color conversions are supported.
return value = LodePNG error code (0 if all went ok, an error if the conversion isn't supported)
The out buffer must have size (w * h * bpp + 7) / 8, where bpp is the bits per pixel
of the output color type (lodepng_get_bpp).
For < 8 bpp images, there should not be padding bits at the end of scanlines.
For 16-bit per channel colors, uses big endian format like PNG does.
Return value is LodePNG error code
*/
unsigned lodepng_convert(unsigned char* out, const unsigned char* in,
	const LodePNGColorMode* mode_out, const LodePNGColorMode* mode_in,
	unsigned w, unsigned h);

#ifdef LODEPNG_COMPILE_DECODER
/*
Settings for the decoder. This contains settings for the PNG and the Zlib
decoder, but not the Info settings from the Info structs.
*/
typedef struct LodePNGDecoderSettings {
	LodePNGDecompressSettings zlibsettings; /*in here is the setting to ignore Adler32 checksums*/

	/* Check LodePNGDecompressSettings for more ignorable errors such as ignore_adler32 */
	unsigned ignore_crc; /*ignore CRC checksums*/
	unsigned ignore_critical; /*ignore unknown critical chunks*/
	unsigned ignore_end; /*ignore issues at end of file if possible (missing IEND chunk, too large chunk, ...)*/
	/* TODO: make a system involving warnings with levels and a strict mode instead. Other potentially recoverable
	   errors: srgb rendering intent value, size of content of ancillary chunks, more than 79 characters for some
	   strings, placement/combination rules for ancillary chunks, crc of unknown chunks, allowed characters
	   in string keys, etc... */

	unsigned color_convert; /*whether to convert the PNG to the color type you want. Default: yes*/

#ifdef LODEPNG_COMPILE_ANCILLARY_CHUNKS
	unsigned read_text_chunks; /*if false but remember_unknown_chunks is true, they're stored in the unknown chunks*/
	/*store all bytes from unknown chunks in the LodePNGInfo (off by default, useful for a png editor)*/
	unsigned remember_unknown_chunks;
#endif /*LODEPNG_COMPILE_ANCILLARY_CHUNKS*/
} LodePNGDecoderSettings;

void lodepng_decoder_settings_init(LodePNGDecoderSettings* settings);
#endif /*LODEPNG_COMPILE_DECODER*/

#ifdef LODEPNG_COMPILE_ENCODER
/*automatically use color type with less bits per pixel if losslessly possible. Default: AUTO*/
typedef enum LodePNGFilterStrategy {
	/*every filter at zero*/
	LFS_ZERO = 0,
	/*every filter at 1, 2, 3 or 4 (paeth), unlike LFS_ZERO not a good choice, but for testing*/
	LFS_ONE = 1,
	LFS_TWO = 2,
	LFS_THREE = 3,
	LFS_FOUR = 4,
	/*Use filter that gives minimum sum, as described in the official PNG filter heuristic.*/
	LFS_MINSUM,
	/*Use the filter type that gives smallest Shannon entropy for this scanline. Depending
	on the image, this is better or worse than minsum.*/
	LFS_ENTROPY,
	/*
	Brute-force-search PNG filters by compressing each filter for each scanline.
	Experimental, very slow, and only rarely gives better compression than MINSUM.
	*/
	LFS_BRUTE_FORCE,
	/*use predefined_filters buffer: you specify the filter type for each scanline*/
	LFS_PREDEFINED
} LodePNGFilterStrategy;

/*Gives characteristics about the integer RGBA colors of the image (count, alpha channel usage, bit depth, ...),
which helps decide which color model to use for encoding.
Used internally by default if "auto_convert" is enabled. Public because it's useful for custom algorithms.*/
typedef struct LodePNGColorStats {
	unsigned colored; /*not grayscale*/
	unsigned key; /*image is not opaque and color key is possible instead of full alpha*/
	unsigned short key_r; /*key values, always as 16-bit, in 8-bit case the byte is duplicated, e.g. 65535 means 255*/
	unsigned short key_g;
	unsigned short key_b;
	unsigned alpha; /*image is not opaque and alpha channel or alpha palette required*/
	unsigned numcolors; /*amount of colors, up to 257. Not valid if bits == 16 or allow_palette is disabled.*/
	unsigned char palette[1024]; /*Remembers up to the first 256 RGBA colors, in no particular order, only valid when numcolors is valid*/
	unsigned bits; /*bits per channel (not for palette). 1,2 or 4 for grayscale only. 16 if 16-bit per channel required.*/
	size_t numpixels;

	/*user settings for computing/using the stats*/
	unsigned allow_palette; /*default 1. if 0, disallow choosing palette colortype in auto_choose_color, and don't count numcolors*/
	unsigned allow_greyscale; /*default 1. if 0, choose RGB or RGBA even if the image only has gray colors*/
} LodePNGColorStats;

void lodepng_color_stats_init(LodePNGColorStats* stats);

/*Get a LodePNGColorStats of the image. The stats must already have been inited.
Returns error code (e.g. alloc fail) or 0 if ok.*/
unsigned lodepng_compute_color_stats(LodePNGColorStats* stats,
	const unsigned char* image, unsigned w, unsigned h,
	const LodePNGColorMode* mode_in);

/*Settings for the encoder.*/
typedef struct LodePNGEncoderSettings {
	LodePNGCompressSettings zlibsettings; /*settings for the zlib encoder, such as window size, ...*/

	unsigned auto_convert; /*automatically choose output PNG color type. Default: true*/

	/*If true, follows the official PNG heuristic: if the PNG uses a palette or lower than
	8 bit depth, set all filters to zero. Otherwise use the filter_strategy. Note that to
	completely follow the official PNG heuristic, filter_palette_zero must be true and
	filter_strategy must be LFS_MINSUM*/
	unsigned filter_palette_zero;
	/*Which filter strategy to use when not using zeroes due to filter_palette_zero.
	Set filter_palette_zero to 0 to ensure always using your chosen strategy. Default: LFS_MINSUM*/
	LodePNGFilterStrategy filter_strategy;
	/*used if filter_strategy is LFS_PREDEFINED. In that case, this must point to a buffer with
	the same length as the amount of scanlines in the image, and each value must <= 5. You
	have to cleanup this buffer, LodePNG will never free it. Don't forget that filter_palette_zero
	must be set to 0 to ensure this is also used on palette or low bitdepth images.*/
	const unsigned char* predefined_filters;

	/*force creating a PLTE chunk if colortype is 2 or 6 (= a suggested palette).
	If colortype is 3, PLTE is _always_ created.*/
	unsigned force_palette;
#ifdef LODEPNG_COMPILE_ANCILLARY_CHUNKS
	/*add LodePNG identifier and version as a text chunk, for debugging*/
	unsigned add_id;
	/*encode text chunks as zTXt chunks instead of tEXt chunks, and use compression in iTXt chunks*/
	unsigned text_compression;
#endif /*LODEPNG_COMPILE_ANCILLARY_CHUNKS*/
} LodePNGEncoderSettings;

void lodepng_encoder_settings_init(LodePNGEncoderSettings* settings);
#endif /*LODEPNG_COMPILE_ENCODER*/


#if defined(LODEPNG_COMPILE_DECODER) || defined(LODEPNG_COMPILE_ENCODER)
/*The settings, state and information for extended encoding and decoding.*/
typedef struct LodePNGState {
#ifdef LODEPNG_COMPILE_DECODER
	LodePNGDecoderSettings decoder; /*the decoding settings*/
#endif /*LODEPNG_COMPILE_DECODER*/
#ifdef LODEPNG_COMPILE_ENCODER
	LodePNGEncoderSettings encoder; /*the encoding settings*/
#endif /*LODEPNG_COMPILE_ENCODER*/
	LodePNGColorMode info_raw; /*specifies the format in which you would like to get the raw pixel buffer*/
	LodePNGInfo info_png; /*info of the PNG image obtained after decoding*/
	unsigned error;
} LodePNGState;

/*init, cleanup and copy functions to use with this struct*/
void lodepng_state_init(LodePNGState* state);
void lodepng_state_cleanup(LodePNGState* state);
void lodepng_state_copy(LodePNGState* dest, const LodePNGState* source);
#endif /* defined(LODEPNG_COMPILE_DECODER) || defined(LODEPNG_COMPILE_ENCODER) */

#ifdef LODEPNG_COMPILE_DECODER
/*
Same as lodepng_decode_memory, but uses a LodePNGState to allow custom settings and
getting much more information about the PNG image and color mode.
*/
unsigned lodepng_decode(unsigned char** out, unsigned* w, unsigned* h,
	LodePNGState* state,
	const unsigned char* in, size_t insize);

/*
Read the PNG header, but not the actual data. This returns only the information
that is in the IHDR chunk of the PNG, such as width, height and color type. The
information is placed in the info_png field of the LodePNGState.
*/
unsigned lodepng_inspect(unsigned* w, unsigned* h,
	LodePNGState* state,
	const unsigned char* in, size_t insize);
#endif /*LODEPNG_COMPILE_DECODER*/

/*
Reads one metadata chunk (other than IHDR) of the PNG file and outputs what it
read in the state. Returns error code on failure.
Use lodepng_inspect first with a new state, then e.g. lodepng_chunk_find_const
to find the desired chunk type, and if non null use lodepng_inspect_chunk (with
chunk_pointer - start_of_file as pos).
Supports most metadata chunks from the PNG standard (gAMA, bKGD, tEXt, ...).
Ignores unsupported, unknown, non-metadata or IHDR chunks (without error).
Requirements: &in[pos] must point to start of a chunk, must use regular
lodepng_inspect first since format of most other chunks depends on IHDR, and if
there is a PLTE chunk, that one must be inspected before tRNS or bKGD.
*/
unsigned lodepng_inspect_chunk(LodePNGState* state, size_t pos,
	const unsigned char* in, size_t insize);

#ifdef LODEPNG_COMPILE_ENCODER
/*This function allocates the out buffer with standard malloc and stores the size in *outsize.*/
unsigned lodepng_encode(unsigned char** out, size_t* outsize,
	const unsigned char* image, unsigned w, unsigned h,
	LodePNGState* state);
#endif /*LODEPNG_COMPILE_ENCODER*/

/*
The lodepng_chunk functions are normally not needed, except to traverse the
unknown chunks stored in the LodePNGInfo struct, or add new ones to it.
It also allows traversing the chunks of an encoded PNG file yourself.

The chunk pointer always points to the beginning of the chunk itself, that is
the first byte of the 4 length bytes.

In the PNG file format, chunks have the following format:
-4 bytes length: length of the data of the chunk in bytes (chunk itself is 12 bytes longer)
-4 bytes chunk type (ASCII a-z,A-Z only, see below)
-length bytes of data (may be 0 bytes if length was 0)
-4 bytes of CRC, computed on chunk name + data

The first chunk starts at the 8th byte of the PNG file, the entire rest of the file
exists out of concatenated chunks with the above format.

PNG standard chunk ASCII naming conventions:
-First byte: uppercase = critical, lowercase = ancillary
-Second byte: uppercase = public, lowercase = private
-Third byte: must be uppercase
-Fourth byte: uppercase = unsafe to copy, lowercase = safe to copy
*/

/*
Gets the length of the data of the chunk. Total chunk length has 12 bytes more.
There must be at least 4 bytes to read from. If the result value is too large,
it may be corrupt data.
*/
unsigned lodepng_chunk_length(const unsigned char* chunk);

/*puts the 4-byte type in null terminated string*/
void lodepng_chunk_type(char type[5], const unsigned char* chunk);

/*check if the type is the given type*/
unsigned char lodepng_chunk_type_equals(const unsigned char* chunk, const char* type);

/*0: it's one of the critical chunk types, 1: it's an ancillary chunk (see PNG standard)*/
unsigned char lodepng_chunk_ancillary(const unsigned char* chunk);

/*0: public, 1: private (see PNG standard)*/
unsigned char lodepng_chunk_private(const unsigned char* chunk);

/*0: the chunk is unsafe to copy, 1: the chunk is safe to copy (see PNG standard)*/
unsigned char lodepng_chunk_safetocopy(const unsigned char* chunk);

/*get pointer to the data of the chunk, where the input points to the header of the chunk*/
unsigned char* lodepng_chunk_data(unsigned char* chunk);
const unsigned char* lodepng_chunk_data_const(const unsigned char* chunk);

/*returns 0 if the crc is correct, 1 if it's incorrect (0 for OK as usual!)*/
unsigned lodepng_chunk_check_crc(const unsigned char* chunk);

/*generates the correct CRC from the data and puts it in the last 4 bytes of the chunk*/
void lodepng_chunk_generate_crc(unsigned char* chunk);

/*
Iterate to next chunks, allows iterating through all chunks of the PNG file.
Input must be at the beginning of a chunk (result of a previous lodepng_chunk_next call,
or the 8th byte of a PNG file which always has the first chunk), or alternatively may
point to the first byte of the PNG file (which is not a chunk but the magic header, the
function will then skip over it and return the first real chunk).
Will output pointer to the start of the next chunk, or at or beyond end of the file if there
is no more chunk after this or possibly if the chunk is corrupt.
Start this process at the 8th byte of the PNG file.
In a non-corrupt PNG file, the last chunk should have name "IEND".
*/
unsigned char* lodepng_chunk_next(unsigned char* chunk, unsigned char* end);
const unsigned char* lodepng_chunk_next_const(const unsigned char* chunk, const unsigned char* end);

/*Finds the first chunk with the given type in the range [chunk, end), or returns NULL if not found.*/
unsigned char* lodepng_chunk_find(unsigned char* chunk, unsigned char* end, const char type[5]);
const unsigned char* lodepng_chunk_find_const(const unsigned char* chunk, const unsigned char* end, const char type[5]);

/*
Appends chunk to the data in out. The given chunk should already have its chunk header.
The out variable and outsize are updated to reflect the new reallocated buffer.
Returns error code (0 if it went ok)
*/
unsigned lodepng_chunk_append(unsigned char** out, size_t* outsize, const unsigned char* chunk);

/*
Appends new chunk to out. The chunk to append is given by giving its length, type
and data separately. The type is a 4-letter string.
The out variable and outsize are updated to reflect the new reallocated buffer.
Returne error code (0 if it went ok)
*/
unsigned lodepng_chunk_create(unsigned char** out, size_t* outsize, unsigned length,
	const char* type, const unsigned char* data);


/*Calculate CRC32 of buffer*/
unsigned lodepng_crc32(const unsigned char* buf, size_t len);
#endif /*LODEPNG_COMPILE_PNG*/


#ifdef LODEPNG_COMPILE_ZLIB
/*
This zlib part can be used independently to zlib compress and decompress a
buffer. It cannot be used to create gzip files however, and it only supports the
part of zlib that is required for PNG, it does not support dictionaries.
*/

#ifdef LODEPNG_COMPILE_DECODER
/*Inflate a buffer. Inflate is the decompression step of deflate. Out buffer must be freed after use.*/
unsigned lodepng_inflate(unsigned char** out, size_t* outsize,
	const unsigned char* in, size_t insize,
	const LodePNGDecompressSettings* settings);

/*
Decompresses Zlib data. Reallocates the out buffer and appends the data. The
data must be according to the zlib specification.
Either, *out must be NULL and *outsize must be 0, or, *out must be a valid
buffer and *outsize its size in bytes. out must be freed by user after usage.
*/
unsigned lodepng_zlib_decompress(unsigned char** out, size_t* outsize,
	const unsigned char* in, size_t insize,
	const LodePNGDecompressSettings* settings);
#endif /*LODEPNG_COMPILE_DECODER*/

#ifdef LODEPNG_COMPILE_ENCODER
/*
Compresses data with Zlib. Reallocates the out buffer and appends the data.
Zlib adds a small header and trailer around the deflate data.
The data is output in the format of the zlib specification.
Either, *out must be NULL and *outsize must be 0, or, *out must be a valid
buffer and *outsize its size in bytes. out must be freed by user after usage.
*/
unsigned lodepng_zlib_compress(unsigned char** out, size_t* outsize,
	const unsigned char* in, size_t insize,
	const LodePNGCompressSettings* settings);

/*
Find length-limited Huffman code for given frequencies. This function is in the
public interface only for tests, it's used internally by lodepng_deflate.
*/
unsigned lodepng_huffman_code_lengths(unsigned* lengths, const unsigned* frequencies,
	size_t numcodes, unsigned maxbitlen);

/*Compress a buffer with deflate. See RFC 1951. Out buffer must be freed after use.*/
unsigned lodepng_deflate(unsigned char** out, size_t* outsize,
	const unsigned char* in, size_t insize,
	const LodePNGCompressSettings* settings);

#endif /*LODEPNG_COMPILE_ENCODER*/
#endif /*LODEPNG_COMPILE_ZLIB*/

#ifdef LODEPNG_COMPILE_DISK
/*
Load a file from disk into buffer. The function allocates the out buffer, and
after usage you should free it.
out: output parameter, contains pointer to loaded buffer.
outsize: output parameter, size of the allocated out buffer
filename: the path to the file to load
return value: error code (0 means ok)
*/
unsigned lodepng_load_file(unsigned char** out, size_t* outsize, const char* filename);

/*
Save a file from buffer to disk. Warning, if it exists, this function overwrites
the file without warning!
buffer: the buffer to write
buffersize: size of the buffer to write
filename: the path to the file to save to
return value: error code (0 means ok)
*/
unsigned lodepng_save_file(const unsigned char* buffer, size_t buffersize, const char* filename);
#endif /*LODEPNG_COMPILE_DISK*/

#ifdef LODEPNG_COMPILE_CPP
/* The LodePNG C++ wrapper uses std::vectors instead of manually allocated memory buffers. */
namespace lodepng {
#ifdef LODEPNG_COMPILE_PNG
	class State: public LodePNGState {
	public:
		State();
		State(const State& other);
		~State();
		State& operator=(const State& other);
	};

#ifdef LODEPNG_COMPILE_DECODER
	/* Same as other lodepng::decode, but using a State for more settings and information. */
	unsigned decode(std::vector<unsigned char>& out, unsigned& w, unsigned& h,
		State& state,
		const unsigned char* in, size_t insize);
	unsigned decode(std::vector<unsigned char>& out, unsigned& w, unsigned& h,
		State& state,
		const std::vector<unsigned char>& in);
#endif /*LODEPNG_COMPILE_DECODER*/

#ifdef LODEPNG_COMPILE_ENCODER
	/* Same as other lodepng::encode, but using a State for more settings and information. */
	unsigned encode(std::vector<unsigned char>& out,
		const unsigned char* in, unsigned w, unsigned h,
		State& state);
	unsigned encode(std::vector<unsigned char>& out,
		const std::vector<unsigned char>& in, unsigned w, unsigned h,
		State& state);
#endif /*LODEPNG_COMPILE_ENCODER*/

#ifdef LODEPNG_COMPILE_DISK
	/*
	Load a file from disk into an std::vector.
	return value: error code (0 means ok)
	*/
	unsigned load_file(std::vector<unsigned char>& buffer, const std::string& filename);

	/*
	Save the binary data in an std::vector to a file on disk. The file is overwritten
	without warning.
	*/
	unsigned save_file(const std::vector<unsigned char>& buffer, const std::string& filename);
#endif /* LODEPNG_COMPILE_DISK */
#endif /* LODEPNG_COMPILE_PNG */

#ifdef LODEPNG_COMPILE_ZLIB
#ifdef LODEPNG_COMPILE_DECODER
	/* Zlib-decompress an unsigned char buffer */
	unsigned decompress(std::vector<unsigned char>& out, const unsigned char* in, size_t insize,
		const LodePNGDecompressSettings& settings = lodepng_default_decompress_settings);

	/* Zlib-decompress an std::vector */
	unsigned decompress(std::vector<unsigned char>& out, const std::vector<unsigned char>& in,
		const LodePNGDecompressSettings& settings = lodepng_default_decompress_settings);
#endif /* LODEPNG_COMPILE_DECODER */

#ifdef LODEPNG_COMPILE_ENCODER
	/* Zlib-compress an unsigned char buffer */
	unsigned compress(std::vector<unsigned char>& out, const unsigned char* in, size_t insize,
		const LodePNGCompressSettings& settings = lodepng_default_compress_settings);

	/* Zlib-compress an std::vector */
	unsigned compress(std::vector<unsigned char>& out, const std::vector<unsigned char>& in,
		const LodePNGCompressSettings& settings = lodepng_default_compress_settings);
#endif /* LODEPNG_COMPILE_ENCODER */
#endif /* LODEPNG_COMPILE_ZLIB */
} /* namespace lodepng */
#endif /*LODEPNG_COMPILE_CPP*/
#endif /*LODEPNG_H inclusion guard*/
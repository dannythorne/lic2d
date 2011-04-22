//##############################################################################
//
// bmp_class.h
//
#ifndef BMP_CLASS_H
#define BMP_CLASS_H

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "bing.h"

#define ROUND floor

#define COUT_RGBS 0

#define SWAP_BYTE_ORDER 0

//#if SWAP_BYTE_ORDER || OSTYPE==darwin
#if SWAP_BYTE_ORDER
// Swap byte order.
#define ENDIAN2(w) ((((w)&0x00ff)<<8)|(((w)&0xff00)>>8))
#define ENDIAN4(w) ((((w)&0x000000ff)<<24)|(((w)&0xff000000)>>24)|(((w)&0x0000ff00)<<8)|(((w)&0x00ff0000)>>8))
#else /* !( SWAP_BYTE_ORDER) */
#define ENDIAN2(w) (w)
#define ENDIAN4(w) (w)
#endif /* SWAP_BYTE_ORDER */

struct bitmap_file_header
{
  // 1 2 bfType 19778 must always be set to 'BM' to declare that this is
  // a .bmp-file.
  char bfType[2];
  // 3 4 bfSize ?? specifies the size of the file in bytes.
  char bfSize[4];
  // 7 2 bfReserved1 0 must always be set to zero.
  char bfReserved1[2];
  // 9 2 bfReserved2 0 must always be set to zero.
  char bfReserved2[2];
  // 11 4 bfOffBits 1078 specifies the offset from the beginning of the
  // file to the bitmap data.
  char bfOffBits[4];

}; /* struct bitmap_file_header */

struct bitmap_info_header
{
  // 15 4 biSize 40 specifies the size of the BITMAPINFOHEADER structure,
  // in bytes.
  char biSize[4];
  // 19 4 biWidth 100 specifies the width of the image, in pixels.
  char biWidth[4];
  // 23 4 biHeight 100 specifies the height of the image, in pixels.
  char biHeight[4];
  // 27 2 biPlanes 1 specifies the number of planes of the target device,
  // must be set to zero. [DT: Should be set to one, right? Not zero.]
  char biPlanes[2];
  // 29 2 biBitCount 8 specifies the number of bits per pixel.
  char biBitCount[2];
  // 31 4 biCompression 0 Specifies the type of compression, usually set
  // to zero (no compression).
  char biCompression[4];
  // 35 4 biSizeImage 0 specifies the size of the image data, in bytes.
  // If there is no compression, it is valid to set this member to zero.
  char biSizeImage[4];
  // 39 4 biXPelsPerMeter 0 specifies the the horizontal pixels per meter
  // on the designated targer device, usually set to zero.
  char biXPelsPerMeter[4];
  // 43 4 biYPelsPerMeter 0 specifies the the vertical pixels per meter
  // on the designated targer device, usually set to zero.
  char biYPelsPerMeter[4];
  // 47 4 biClrUsed 0 specifies the number of colors used in the bitmap,
  // if set to zero the number of colors is calculated using the biBitCount
  // member.
  char biClrUsed[4];
  // 51 4 biClrImportant 0 specifies the number of color that are
  // 'important' for the bitmap, if set to zero, all colors are important.
  char biClrImportant[4];

}; /* struct bitmap_info_header */

struct bmp_hdr_struct
{
  char type[2];  // 'BM'
  char size[4];  // Size of file in bytes.
  char rsvd[4];  // 0
  char dptr[4];  // Offset of data in bytes.

  char forty[4]; // Size of info header (=40).
  char width[4]; // Width of image in pixels
  char height[4];// Height of image in pixels.
  char planes[2];// Must be set to one.
  char depth[2]; // Color depth: bits per pixel.
  char zeros[24];// Unused parameters.

}; /* struct bmp_hdr_struct */

struct rgb_quad
{
  // 1 1 rgbBlue - specifies the blue part of the color.
  char Blue;
  // 2 1 rgbGreen - specifies the green part of the color.
  char Green;
  // 3 1 rgbRed - specifies the red part of the color.
  char Red;
  // 4 1 rgbReserved - must always be set to zero.
  char Reserved;
};

class bmp_class
{
private:

  std::string filename;
  std::ifstream fin;

  char type[2];  // 'BM'
  int size;  // Size of file in bytes.
  int rsvd;  // 0
  int dptr;  // Offset of data in bytes.

  int forty; // Size of info header (=40).
  int width; // Width of image in pixels
  int height;// Height of image in pixels.
  short int planes;// Must be set to one.
  short int depth; // Color depth: bits per pixel.
  int zeros[6];// Unused parameters.

  int *r;
  int *g;
  int *b;

public:

  bmp_class()
  {
    r = NULL;
    g = NULL;
    b = NULL;
  }

  bmp_class( std::string arg_fn)
  {
    filename = arg_fn;

    std::cout << "Opening \"" << filename << "\"" << std::endl;
    fin.open( filename.data());

    read_header();
    disp_header();

    r = new int[ width*height];
    g = new int[ width*height];
    b = new int[ width*height];

    switch( depth)
    {
      case 1:
        BING("Mono.");
        //read_mono_color();
        break;
      case 4:
        BING("Four.");
        //read_four_color();
        break;
      case 8:
        BING("What.");
        //read_what_color();
        break;
      case 24:
        BING("True.");
        read_true_color();
        break;
      default:
        BING("Unhandled case.");
        break;
    }

  } /* bmp_class( std::string arg_fn) */

  ~bmp_class()
  {
    if( filename.length() > 0)
    {
      std::cout << "Closing \"" << filename << "\"" << std::endl;
    }
    if( fin.is_open())
    {
      fin.close();
    }

    if( r) { delete [] r;}
    if( g) { delete [] g;}
    if( b) { delete [] b;}

  }

public:

  void read_header();
  void disp_header();

  void read_true_color();

  void create(
    const std::string arg_filename,
    const int arg_width,
    const int arg_height,
    const int arg_depth = 24 );

  void init_header();

  void write();

  void set_xyrgb(
    const int x,
    const int y,
    const int arg_r,
    const int arg_g,
    const int arg_b );

  void set_filename( const std::string arg_filename)
  {
    filename = arg_filename;
  }

  void flipud();

  void scroll_left( const int num_cols);
  void scroll_down( const int num_rows);

public: // Accessors.

  // Bytes of pixel information per row of the bitmap.
  int bytes_of_pixels_per_row() const
  {
    return (int)ceil( ( ((double)width)*( (double)depth))/8.);
  }

  // Bitmaps pad rows to preserve 4-byte boundaries.
  // The length of a row in the file will be
  //
  //   bytes_of_pixels_per_row + bytes_of_pad .
  //
  int bytes_of_pad() const
  {
    return ( 4 - bytes_of_pixels_per_row() % 4) % 4;
  }

  int num_pixels() const
  {
    return width*height;
  }

  int get_width() const { return width;}
  int get_height() const { return height;}

  int get_r( const int x, const int y) const { return r[x+y*width];}
  int get_g( const int x, const int y) const { return g[x+y*width];}
  int get_b( const int x, const int y) const { return b[x+y*width];}

}; /* class bmp_class */

#endif

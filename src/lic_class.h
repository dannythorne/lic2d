//##############################################################################
//
// lic_class.h
//

#include <iostream>
#include <cstdlib>

#include "bmp_class.h"

#define PI 3.141592654
#define EPS 1e-5

namespace dthorne0_lic_class_ns
{
class lic_class
{
public:

  typedef double real;
  typedef float value_type;
  typedef unsigned char texture_type;

private: // Member variables

  int ni;
  int nj;
  int ns; // Maximum number of nodes in streamline is 2*ns+1.

  value_type ***v; // Raw vector data.
  texture_type **r; // Input texture, red channel.
  texture_type **g; // Input texture, green channel.
  texture_type **b; // Input texture, blue channel.
  texture_type **rprime; // Input texture.
  texture_type **gprime; // Input texture.
  texture_type **bprime; // Input texture.

  value_type **sfore; // Stream line coordinates.
  value_type **sback; // Stream line coordinates.
  value_type *dfore; // Stream line length.
  value_type *dback; // Stream line length.
  int **tfore; // Stream pixel coordinates.
  int **tback; // Stream pixel coordinates.

  texture_type texture_style_index;

  real vxmin, vxmax;
  real vymin, vymax;
  real vmin, vmax;

  bool periodic;

public: // Constructor(s)

  lic_class(
    texture_type arg_texture_style_index = 3,
    const char *filename = NULL,
    const int arg_ns = 16);
  ~lic_class();

public: // Accessors

  void display( std::ostream &o) const;

public: // Modifiers

  void set_texture_style_index( const texture_type n) { texture_style_index=n;}

  void compute_fprime_approx();
  void compute_fprime();
  void compute_streamline_pixels(
         int i0, int j0,
         int &countpos, int &countneg );
  void compute_streamline(
         const int i0, const int j0,
         int &countpos, int &countneg );
  value_type compute_h_fore_approx( const int n0, const int countpos);
  value_type compute_h_back_approx( const int n0, const int countneg);
  value_type compute_h_fore( const int n);
  value_type compute_h_back( const int n);

private: // Helpers
  bool are_valid_indices( const int i, const int j) const;
  int which_quadrant_vector_points_into(
        const real vx, const real vy) const;
  bool on_bottom_wall( const real x, const real y) const { /*Don't need.*/}
  bool on_left_wall( const real x, const real y) const { /*Don't need.*/}
  bool pointing_at_top_wall(
    const real x,
    const real y,
    const real vx,
    const real vy ) const;
  bool pointing_at_bottom_wall(
    const real x,
    const real y,
    const real vx,
    const real vy ) const;
  bool pointing_at_right_wall(
    const real x,
    const real y,
    const real vx,
    const real vy ) const;
  bool pointing_at_left_wall(
    const real x,
    const real y,
    const real vx,
    const real vy ) const;
  bool not_off_domain( real &x, real &y) const;

  texture_type texture_noise( const int i, const int j)
  { // Standard white noise: P(black)==P(white)==.5
    return ( real(rand())/real(RAND_MAX) < .5)?(0):(255);
  }
  void texture_velocity_blue_green_white( const int i, const int j)
  {
    real u;
    texture_type ut;

    u = real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                  + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)/real(vmax-vmin);

    ut = (texture_type)floor( 255*u);

    g[j][i] = (i%4||!(j%4))
            ?((j%4||!(i%4))?(0):(3*int(ut)/4+rand()%(int(ut)/4+1)))
            :(3*int(ut)/4+rand()%(int(ut)/4+1));

    r[j][i] = 0;

    ut = (texture_type)floor( 255.-255*u);

    b[j][i] = (i%4||!(j%4))
            ?((j%4||!(i%4))?(0):(3*int(ut)/4+rand()%(int(ut)/4+1)))
            :(3*int(ut)/4+rand()%(int(ut)/4+1));
    if( 0
      ||(!((i+1)%4) && !((j+1)%4))
      ||(!((i+1)%4) && !((j+2)%4))
      ||(!((i+1)%4) && !((j+3)%4))
      ||(!((i+2)%4) && !((j+1)%4))
      ||(!((i+2)%4) && !((j+2)%4))
      ||(!((i+2)%4) && !((j+3)%4))
      ||(!((i+3)%4) && !((j+1)%4))
      ||(!((i+3)%4) && !((j+2)%4))
      ||(!((i+3)%4) && !((j+3)%4))
      )
    {
        ut = (texture_type)floor( 223*u);
        r[j][i] = 32+rand()%(int(ut)+1);
        g[j][i] = r[j][i]-32 + rand()%33;
        b[j][i] = r[j][i]-32 + rand()%33;
    }

  }

  void texture_velocity_red_green_noise( const int i, const int j)
  {
    real u;
    texture_type ut;

    u = real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                  + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)/real(vmax-vmin);
    ut = (texture_type)floor(255*u);
    b[j][i] = (real(rand())/real(RAND_MAX)<.1)?(0):(rand()%(int(ut)+1));

    u = (real( fabs(v[j][i][0]))-vxmin)/real(vxmax-vxmin);
    ut = (texture_type)floor(255*u);
    g[j][i] = (real(rand())/real(RAND_MAX)<.1)?(0):(rand()%(int(ut)+1));
    //g[j][i] = r[j][i];

    u = (real( fabs(v[j][i][1]))-vymin)/real(vymax-vymin);
    ut = (texture_type)floor(255*u);
    r[j][i] = (real(rand())/real(RAND_MAX)<.1)?(0):(rand()%(int(ut)+1));
  }

  void texture_red_green_blue_noise( const int i, const int j)
  {
      r[j][i] = (i%4)?(0):(128+rand()%128);
      g[j][i] = (j%4||!(i%4))?(0):(128+rand()%128);
      b[j][i] = ( (!((i+1)%4) && !((j+1)%4))
                ||(!((i+1)%4) && !((j+3)%4))
                ||(!((i+2)%4) && !((j+2)%4))
                ||(!((i+3)%4) && !((j+1)%4))
                ||(!((i+3)%4) && !((j+3)%4)) )?(128+rand()%128):(0);
      if( (!((i+2)%4) && !((j+1)%4))
        ||(!((i+1)%4) && !((j+2)%4))
        ||(!((i+2)%4) && !((j+3)%4))
        ||(!((i+3)%4) && !((j+2)%4)) )
      {
        r[j][i] = rand()%128;
        g[j][i] = r[j][i];
        b[j][i] = r[j][i];
      }
  }
  void texture_blue_green( const int i, const int j)
  {
    real u;
    texture_type ut;

    r[j][i] = 0;
#if 0
    u = real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                  + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)/real(vmax-vmin);
    ut = (texture_type)floor(255*u);
    r[j][i] = (real(rand())/real(RAND_MAX)<.1)?(0):(rand()%(int(ut)+1));
#endif

    u = (real( fabs(v[j][i][0]))-vxmin)/real(vxmax-vxmin);
    ut = (texture_type)floor(255*u);
    g[j][i] = (real(rand())/real(RAND_MAX)<.1)?(0):(rand()%(int(ut)+1));

    u = (real( fabs(v[j][i][1]))-vymin)/real(vymax-vymin);
    ut = (texture_type)floor(255*u);
    b[j][i] = (real(rand())/real(RAND_MAX)<.1)?(0):(rand()%(int(ut)+1));
  }

}; //class lic_class

std::ostream &operator<<( std::ostream &o, const lic_class &arg_lic);

} //namespace dthorne0_lic_class_ns

// vim: foldmethod=syntax foldlevel=1

//##############################################################################
//
// lic_class.cc
//

#include "lic_class.h"
#include "bmp_class.h"

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#define DETECT_LOOPS 1

namespace dthorne0_lic_class_ns
{
lic_class::lic_class(
  texture_type arg_texture_style_index,
  const char *filename,
  const int arg_ns)
{
  int i, j;
  real u;

  ni = 512;
  nj = 512;
  ns = arg_ns;

  vxmin=999999.; vxmax=0.;
  vymin=999999.; vymax=0.;
  vmin=999999.; vmax=0.;

  texture_style_index = arg_texture_style_index;
  periodic = true;

  // Initialize velocity field.
  if( filename!=NULL)
  {
    ifstream in;
    in.open( filename);
    cout << "\"" << filename << "\" successfully opened." << endl;
    in >> ni;
    in >> nj;
    v = new value_type**[nj];
    for( j=0; j<nj; j++)
    {
      v[j] = new value_type*[ni];
      for( i=0; i<ni; i++)
      {
        v[j][i] = new value_type[2];
        in >> v[j][i][0];
        if( vxmin > fabs(v[j][i][0])) { vxmin = fabs(v[j][i][0]);}
        if( vxmax < fabs(v[j][i][0])) { vxmax = fabs(v[j][i][0]);}
        in >> v[j][i][1];
        if( vymin > fabs(v[j][i][1])) { vymin = fabs(v[j][i][1]);}
        if( vymax < fabs(v[j][i][1])) { vymax = fabs(v[j][i][1]);}

        u = sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                + fabs(v[j][i][1])*fabs(v[j][i][1]));
        if( vmin > u) { vmin = u;}
        if( vmax < u) { vmax = u;}

      }
    }
    in.close();
  }
  else
  {
    // Initialize vector data.
    v = new value_type**[nj];
    value_type squareroot = 0.0;
    value_type x, y;
    value_type theta;
    cout << "Default circular velocity field: "
         << ni << "x" << nj << "." << endl;
    for( j=0; j<nj; j++)
    {
      v[j] = new value_type*[ni];

      y = j - ( value_type(nj)/2. - 1./2.);

      for( i=0; i<ni; i++)
      {
        v[j][i] = new value_type[2];

        // Initial values.
        x = i - ( value_type(ni)/2. - 1./2.);
        squareroot = sqrt(x*x + y*y);
        if( squareroot != 0)
        {
          // Flow in a circle around center of domain.
          switch( which_quadrant_vector_points_into(x,y))
          {
            case 1: theta = atan(y/x); break;
            case 2: theta = PI+atan(y/x); break;
            case 3: theta = PI+atan(y/x); break;
            case 4: theta = 2*PI+atan(y/x); break;
            default:
              cout << "ERROR: Unhandled case: quadrant "
                   << which_quadrant_vector_points_into(x,y) << endl;
              exit(1);
          }
          v[j][i][0] = /*sin(.5*(theta+PI/3.))**/(value_type)(-y)/squareroot;
          v[j][i][1] = /*sin(.5*(theta+PI/3.))**/(value_type)x/squareroot;
        }
        else
        {
          v[j][i][0] = 0.;
          v[j][i][1] = 0.;
        }
        if( vxmin > fabs(v[j][i][0])) { vxmin = fabs(v[j][i][0]);}
        if( vxmax < fabs(v[j][i][0])) { vxmax = fabs(v[j][i][0]);}
        if( vymin > fabs(v[j][i][1])) { vymin = fabs(v[j][i][1]);}
        if( vymax < fabs(v[j][i][1])) { vymax = fabs(v[j][i][1]);}

        u = sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                + fabs(v[j][i][1])*fabs(v[j][i][1]));
        if( vmin > u) { vmin = u;}
        if( vmax < u) { vmax = u;}

#if 1 // Flow uniformly to the right.
        if( i==0 && j==0)
        {
    cout << "Overriding circular velocity field with uniform flow: "
         << ni << "x" << nj << "." << endl;
        }
          v[j][i][0] = 0.;
          v[j][i][1] = 1.;
#endif
      }
    }
  }

  // Initialize input texture.
  cout << "Input texture #"
       << (int)texture_style_index << "." << endl;
  r = new texture_type*[nj];
  g = new texture_type*[nj];
  b = new texture_type*[nj];
  for( j=0; j<nj; j++)
  {
    r[j] = new texture_type[ni];
    g[j] = new texture_type[ni];
    b[j] = new texture_type[ni];
    for( i=0; i<ni; i++)
    {
      switch( (int)texture_style_index)
      {
        case 1: { // Black and white noise.
          r[j][i] = texture_noise( i, j);
          g[j][i] = r[j][i];
          b[j][i] = r[j][i];
          break;
        }
        case 2: { // Color noise.
          r[j][i] = texture_noise( i, j);
          g[j][i] = texture_noise( i, j);
          b[j][i] = texture_noise( i, j);
          break;
        }
        case 3: { // Velocity based. Blue/green frame with white-noise fill.
          texture_velocity_blue_green_white( i, j);
          break;
        }
        case 4: { // Velocity based. Red/green noise.
          texture_velocity_red_green_noise( i, j);
          break;
        }
        case 5: { // Red/Green noise frame with blue noise fill.
          texture_red_green_blue_noise( i, j);
          break;
        }
        case 6: { // Velocity based. Blue/Green
          texture_blue_green( i, j);
          break;
        }
        default:
          cout << "ERROR: Unhandled case." << endl;
          break;
      }

#if 0
#if 0
      //r[j][i] = 128*((i%2)^(j%2));
      //r[j][i] = 128*((i%2)+(1-(j%2)));
      //r[j][i] = rand()%256;

    //r[j][i] = (real(rand())/real(RAND_MAX)<.5)
    //         ?((real(rand())/real(RAND_MAX)<.5)?(0):(255))
    //         :(rand()%256);

      //r[j][i] = (i%4)?(0):(255);
      //g[j][i] = (j%4||!(i%4))?(0):(255);
      //b[j][i] = ( (!((i+1)%4) && !((j+1)%4))
      //          ||(!((i+1)%4) && !((j+3)%4))
      //          ||(!((i+2)%4) && !((j+2)%4))
      //          ||(!((i+3)%4) && !((j+1)%4))
      //          ||(!((i+3)%4) && !((j+3)%4)) )?(255):(0);
#else
#if 0
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
#else
      //g[j][i] = (i%4||!(j%4))?(0)
      //:(int(floor(255.*real(fabs(v[j][i][1])-vymin)/real(vymax-vymin))));
      //b[j][i] = (j%4||!(i%4))?(0)
      //:(int(floor(255.*real(fabs(v[j][i][0])-vxmin)/real(vxmax-vxmin))));
      u =
      (   (floor(255.*(0.+real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                                + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)
                     /real(vmax-vmin)))));
      //g[j][i] = (i%4||!(j%4))?( (j%4||!(i%4))?(0):(int(u)/2+rand()%(int(u)/2+1)) ):(int(u)/2+rand()%(int(u)/2+1));
      g[j][i] = (i%4||!(j%4))?( (j%4||!(i%4))?(0):(3*int(u)/4+rand()%(int(u)/4+1)) ):(3*int(u)/4+rand()%(int(u)/4+1));
      //g[j][i] = (i%4||!(j%4))?( (j%4||!(i%4))?(0):(int(u)) ):(int(u));

      //r[j][i] = (i%4||!(j%4))?( (j%4||!(i%4))?(0):(int(u)) ):(int(u));
      r[j][i] = 0;//(i%4||!(j%4))?( (j%4||!(i%4))?(0):(int(u)) ):(int(u));
      //r[j][i] = (u>200)?(3*int(u)/4+rand()%(int(u)/4+1)):(0);

      u =
      (   (floor(255.*(1.-real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                                + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)
                     /real(vmax-vmin)))));
      //b[j][i] = (i%4||!(j%4))?( (j%4||!(i%4))?(0):(int(u)/2+rand()%(int(u)/2+1)) ):(int(u)/2+rand()%(int(u)/2+1));
      b[j][i] = (i%4||!(j%4))?( (j%4||!(i%4))?(0):(3*int(u)/4+rand()%(int(u)/4+1)) ):(3*int(u)/4+rand()%(int(u)/4+1));
      //b[j][i] = (i%4||!(j%4))?( (j%4||!(i%4))?(0):(int(u)) ):(int(u));

      //r[j][i] = ( (!((i+1)%4) && !((j+1)%4))
      //          ||(!((i+1)%4) && !((j+3)%4))
      //          ||(!((i+2)%4) && !((j+2)%4))
      //          ||(!((i+3)%4) && !((j+1)%4))
      //          ||(!((i+3)%4) && !((j+3)%4)) )
      //?(int(floor(255.*real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
      //                           + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)
      //                /real(vmax-vmin))))
      //:(0);
      if( (!((i+2)%4) && !((j+1)%4))
        ||(!((i+1)%4) && !((j+2)%4))
        ||(!((i+2)%4) && !((j+3)%4))
        ||(!((i+3)%4) && !((j+2)%4))
        ||(!((i+1)%4) && !((j+1)%4))
        ||(!((i+1)%4) && !((j+3)%4))
        ||(!((i+2)%4) && !((j+2)%4))
        ||(!((i+3)%4) && !((j+1)%4))
        ||(!((i+3)%4) && !((j+3)%4)) )
      {
#if 1
#if 1
        u =
        (   (floor(223.*(0.+real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                                  + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)
                       /real(vmax-vmin)))));
        r[j][i] = 32+rand()%(int(u)+1);
        g[j][i] = r[j][i]-32 + rand()%33;
        b[j][i] = r[j][i]-32 + rand()%33;
#else
        u =
        (   (floor(255.*(0.+real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                                  + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)
                       /real(vmax-vmin)))));
        r[j][i] = rand()%(int(u)+1);
        u =
        (   (floor(255.*(1.-real( sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
                                  + fabs(v[j][i][1])*fabs(v[j][i][1]))-vmin)
                       /real(vmax-vmin)))));
        b[j][i] = rand()%(int(u)+1);
        g[j][i] = 0;
#endif
#else
        r[j][i] = 0;
        g[j][i] = 0;
        b[j][i] = 0;
#endif
      }
#endif
#endif
#endif

      //r[j][i] = (i%2)?(0):(rand()%256);
      //g[j][i] = (j%2)?(0):(rand()%256);
      //b[j][i] = ((i+1)%2 || (j+1)%2)?(0):(rand()%256);

      //r[j][i] = (int)floor(255.*(real(i)/real(ni)));
      //r[j][i] = (int)floor(255.*(real(j)/real(nj)));
    }
  }

  // Initialize output texture.
  rprime = new texture_type*[nj];
  gprime = new texture_type*[nj];
  bprime = new texture_type*[nj];
  for( j=0; j<nj; j++)
  {
    rprime[j] = new texture_type[ni];
    gprime[j] = new texture_type[ni];
    bprime[j] = new texture_type[ni];
    for( i=0; i<ni; i++)
    {
      rprime[j][i] = 0;
      gprime[j][i] = 0;
      bprime[j][i] = 0;
    }
  }

  // Initialize stream line pixels.
  tfore = new int*[ns+1];
  tback = new int*[ns+1];
  for( i=0; i<ns+1; i++)
  {
    tfore[i] = new int[2];
    tback[i] = new int[2];

    tfore[i][0] = 0;
    tback[i][0] = 0;

    tfore[i][1] = 0;
    tback[i][1] = 0;
  }

  // Initialize stream line coordinates.
  sfore = new value_type*[ns+1];
  sback = new value_type*[ns+1];
  for( i=0; i<ns+1; i++)
  {
    sfore[i] = new value_type[2];
    sback[i] = new value_type[2];

    sfore[i][0] = 0;
    sback[i][0] = 0;

    sfore[i][1] = 0;
    sback[i][1] = 0;
  }

  // Initialize stream line length.
  dfore = new value_type[ns+1];
  dback = new value_type[ns+1];

  //compute_fprime_approx();
#if 0
  int countpos, countneg;
  int seedi=ni/4;
  int seedj=nj/4;
#if 0
  for( seedj=0; seedj<nj-1; seedj++)
  {
    for( seedi=0; seedi<ni-1; seedi++)
    {
      compute_streamline( seedi, seedj, countpos, countneg);
      for( i=0; i<countpos; i++) { rprime[ tfore[i][1] ][ tfore[i][0] ] = i;}
      for( i=0; i<countneg; i++) { rprime[ tback[i][1] ][ tback[i][0] ] = i;}
      rprime[ seedj][ seedi] = 0;
    }
  }
#else
  compute_streamline( seedi, seedj, countpos, countneg);
#endif
  cout << "countpos = " << countpos << endl;
  cout << "countneg = " << countneg << endl;
#if 1
  cout << "sfore: ";
  for( i=0; i<=countpos; i++)
  {
    cout << sfore[i][0] << " " << sfore[i][1] << " -- ";
  }
  cout << endl;
  cout << "tfore: ";
  for( i=0; i<countpos; i++)
  {
    cout << tfore[i][0] << " " << tfore[i][1] << " -- ";
    rprime[ tfore[i][1] ][ tfore[i][0] ] = 128;
  }
  cout << endl;
#endif
#if 1
  cout << "sback: ";
  for( i=0; i<=countneg; i++)
  {
    cout << sback[i][0] << " " << sback[i][1] << " -- ";
  }
  cout << endl;
  cout << "tback: ";
  for( i=0; i<countneg; i++)
  {
    cout << tback[i][0] << " " << tback[i][1] << " -- ";
    rprime[ tback[i][1] ][ tback[i][0] ] = 255;
  }
  cout << endl;
#endif
  rprime[ seedj][ seedi] = 64;
#else
  compute_fprime();
#endif

} // lic_class::lic_class()

void lic_class::display( std::ostream &o) const
{
  int i, j;

#if 0
  o << "Vector data:" << endl;
  for( j=nj-1; j>=0; j--)
  {
    cout << "  ";
    for( i=0; i<ni; i++)
    {
      cout << "("
           << setprecision(1) << setw(4) << v[j][i][0] << ","
           << setw(4)<< v[j][i][1] << ")";
    }
    cout << endl;
  }
#endif

#if 0
  o << "Input texture:" << endl;
  for( j=nj-1; j>=0; j--)
  {
    cout << "  ";
    for( i=0; i<ni; i++)
    {
      cout << setw(3)<< int(r[j][i]) << " ";
    }
    cout << endl;
  }
#endif

#if 0
  o << "Output texture:" << endl;
  for( j=nj-1; j>=0; j--)
  {
    cout << "  ";
    for( i=0; i<ni; i++)
    {
      cout << setw(3)<< int(rprime[j][i]) << " ";
    }
    cout << endl;
  }
#endif

  bmp_class bmp_in;
  bmp_in.create("lic_in.bmp",ni,nj);
  for( j=nj-1; j>=0; j--)
  {
    for( i=0; i<ni; i++)
    {
      bmp_in.set_xyrgb( i, j, (r[j][i]), (g[j][i]), (b[j][i]));
    }
  }
  //bmp_in.scroll_left(128);
  //bmp_in.scroll_down(64);
  bmp_in.write();

  bmp_class bmp;
  bmp.create("lic.bmp",ni,nj);
  int min_val = 255;
  int max_val = 0;
  double u;
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
      if( min_val > rprime[j][i]) { min_val = rprime[j][i];}
      if( max_val < rprime[j][i]) { max_val = rprime[j][i];}

      if( min_val > gprime[j][i]) { min_val = gprime[j][i];}
      if( max_val < gprime[j][i]) { max_val = gprime[j][i];}

      if( min_val > bprime[j][i]) { min_val = bprime[j][i];}
      if( max_val < bprime[j][i]) { max_val = bprime[j][i];}
    }
  }
  //cout << "(min_val,max_val)=(" << min_val << "," << max_val << ")" << endl;
  for( j=nj-1; j>=0; j--)
  {
    for( i=0; i<ni; i++)
    {
      //bmp.set_xyrgb( i, j,
      //  (int)floor(255.*rprime[j][i]),
      //  (int)floor(255.*rprime[j][i]),
      //  (int)floor(255.*rprime[j][i]));
//cout
//<< (int)rprime[j][i]
//<< " --> "
//<< (int)floor(255.*(real(rprime[j][i]-min_val)/real(max_val-min_val)))
//<< endl;

      u = sqrt( fabs(v[j][i][0])*fabs(v[j][i][0])
              + fabs(v[j][i][1])*fabs(v[j][i][1]));

      if( u!=0.)
      {
        bmp.set_xyrgb( i, j,
          (int)floor(255.*(real(rprime[j][i]-min_val)/real(max_val-min_val))),
          (int)floor(255.*(real(gprime[j][i]-min_val)/real(max_val-min_val))),
          (int)floor(255.*(real(bprime[j][i]-min_val)/real(max_val-min_val))) );
      }
      else
      {
        bmp.set_xyrgb( i, j, 255, 255, 255);
      }
    }
  }

#if 0
  cout << __FILE__ << " " << __LINE__ << " >> "
       << "NOTE: Shifting domain for better viewing of edges." << endl;
  bmp.scroll_left(128);
  bmp.scroll_down(64);
#endif
  bmp.write();

} // void lic_class::display( std::ostream &o) const

void lic_class::compute_fprime()
{
  int i, j, n;
  int countpos, countneg;
  value_type h;
  value_type rh_sum_fore, rh_sum_back;
  value_type gh_sum_fore, gh_sum_back;
  value_type bh_sum_fore, bh_sum_back;
  value_type h_sum_fore;
  value_type h_sum_back;

  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
      // Compute stream line through (i+.5,j+.5).
      compute_streamline( i, j, countpos, countneg);

      // Compute convolutions (weights) and sums.
      rh_sum_fore = 0;
      gh_sum_fore = 0;
      bh_sum_fore = 0;
      h_sum_fore = 0;

      n = 0;
      h = compute_h_fore( n);
      rh_sum_fore += r[ tfore[n][1]][ tfore[n][0]]*h;
      gh_sum_fore += g[ tfore[n][1]][ tfore[n][0]]*h;
      bh_sum_fore += b[ tfore[n][1]][ tfore[n][0]]*h;
      h_sum_fore += h;
      for( n=1; n<countpos-1; n++)
      {
        h = compute_h_fore( n);
        rh_sum_fore += r[ tfore[n][1]][ tfore[n][0]]*h;
        gh_sum_fore += g[ tfore[n][1]][ tfore[n][0]]*h;
        bh_sum_fore += b[ tfore[n][1]][ tfore[n][0]]*h;
        h_sum_fore += h;
      }

      rh_sum_back = 0;
      gh_sum_back = 0;
      bh_sum_back = 0;
      h_sum_back = 0;

      n = 0;
      h = compute_h_back( n);
      rh_sum_back += r[ tback[n][1]][ tback[n][0]]*h;
      gh_sum_back += g[ tback[n][1]][ tback[n][0]]*h;
      bh_sum_back += b[ tback[n][1]][ tback[n][0]]*h;
      h_sum_back += h;
      for( n=1; n<countneg-1; n++)
      {
        h = compute_h_back( n);
        rh_sum_back += r[ tback[n][1]][ tback[n][0]]*h;
        gh_sum_back += g[ tback[n][1]][ tback[n][0]]*h;
        bh_sum_back += b[ tback[n][1]][ tback[n][0]]*h;
        h_sum_back += h;
      }

      // Compute fprime(i,j)
      rprime[j][i] = texture_type(
          (int)floor( ( rh_sum_fore
              + rh_sum_back ) / ( h_sum_fore
                + h_sum_back)));
      gprime[j][i] = texture_type(
          (int)floor( ( gh_sum_fore
              + gh_sum_back ) / ( h_sum_fore
                + h_sum_back)));
      bprime[j][i] = texture_type(
          (int)floor( ( bh_sum_fore
              + bh_sum_back ) / ( h_sum_fore
                + h_sum_back)) );

#if 1
      if( countpos < ns/2 && countneg < ns/2)
        //if( countpos <= 1 && countneg <= 1)
      {
        //rprime[j][i] = 0;
        //gprime[j][i] = 0;
        //bprime[j][i] = 0;
        rprime[j][i] = (int)floor(rprime[j][i]*double(countpos+countneg)/double(1*(ns+0)));
        gprime[j][i] = (int)floor(gprime[j][i]*double(countpos+countneg)/double(1*(ns+0)));
        bprime[j][i] = (int)floor(bprime[j][i]*double(countpos+countneg)/double(1*(ns+0)));
      }
#endif

    }
  }

}

void lic_class::compute_fprime_approx()
{
  int i, j, n;
  int countpos, countneg;
  value_type h;
  value_type fh_sum_fore;
  value_type fh_sum_back;
  value_type h_sum_fore;
  value_type h_sum_back;

  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
      // Compute stream line through (i,j)th pixel.
      compute_streamline_pixels( i, j, countpos, countneg);

      // Compute convolutions (weights) and sums.
      fh_sum_fore = 0;
      h_sum_fore = 0;

      n = 0;
      h = .5;
      fh_sum_fore += r[ tfore[n][1]][ tfore[n][0]]*h;
      h_sum_fore += h;
      for( n=1; n<=countpos; n++)
      {
        h = compute_h_fore_approx( n, countpos);
        fh_sum_fore += r[ tfore[n][1]][ tfore[n][0]]*h;
        h_sum_fore += h;
      }

      fh_sum_back = 0;
      h_sum_back = 0;

      n = 0;
      h = .5;
      fh_sum_back += r[ tback[n][1]][ tback[n][0]]*h;
      h_sum_back += h;
      for( n=1; n<=countneg; n++)
      {
        h = compute_h_back_approx( n, countneg);
        fh_sum_back += r[ tback[n][1]][ tback[n][0]]*h;
        h_sum_back += h;
      }

      // Compute rprime(i,j)
      rprime[j][i] = texture_type(
        (int)floor( ( fh_sum_fore
                    + fh_sum_back ) / ( h_sum_fore
                                      + h_sum_back)));

    }
  }

}

void lic_class::compute_streamline_pixels(
       int i0, int j0,
       int &countpos, int &countneg )
{
  int i, j;
  bool still_inside_domain;
  real vx, vy, theta;

  // Compute streamline forward.
  countpos = 0;
  tfore[countpos][0] = i0;
  tfore[countpos][1] = j0;
  i = i0;
  j = j0;
  still_inside_domain = true;
  while( countpos < ns && still_inside_domain)
  {
    vx = v[j][i][0];
    vy = v[j][i][1];

    theta = fabs( atan( vy/vx));
         if( vx >= 0 && vy >= 0) { if( theta < PI/4.) { i++;} else { j++;}}
    else if( vx >= 0 && vy <  0) { if( theta < PI/4.) { i++;} else { j--;}}
    else if( vx <  0 && vy >= 0) { if( theta < PI/4.) { i--;} else { j++;}}
    else if( vx <  0 && vy <  0) { if( theta < PI/4.) { i--;} else { j--;}}

    if( are_valid_indices( i, j))
    {
      countpos++;
      tfore[countpos][0] = i;
      tfore[countpos][1] = j;
    }
    else
    {
      still_inside_domain = false;
    }
  }

  // Compute streamline backward.
  countneg = 0;
  tback[countneg][0] = i0;
  tback[countneg][1] = j0;
  i = i0;
  j = j0;
  still_inside_domain = true;
  while( countneg < ns && still_inside_domain)
  {
    vx = -v[j][i][0];
    vy = -v[j][i][1];

    theta = fabs( atan( vy/vx));
         if( vx >= 0 && vy >= 0) { if( theta < PI/4.) { i++;} else { j++;}}
    else if( vx >= 0 && vy <  0) { if( theta < PI/4.) { i++;} else { j--;}}
    else if( vx <  0 && vy >= 0) { if( theta < PI/4.) { i--;} else { j++;}}
    else if( vx <  0 && vy <  0) { if( theta < PI/4.) { i--;} else { j--;}}

    if( are_valid_indices( i, j))
    {
      countneg++;
      tback[countneg][0] = i;
      tback[countneg][1] = j;
    }
    else
    {
      still_inside_domain = false;
    }
  }

} // void lic_class::compute_streamline_pixels(

void lic_class::compute_streamline(
       const int i0, const int j0,
       int &countpos, int &countneg )
{
  int i;
  int j;
  int n;
  real x, ds_x;
  real y, ds_y;
  real vx;
  real vy;
  int k;
  bool loop;

 if(1)
 {
  //############################################################################
  //
  // F O R W A R D
  //
  i = i0;
  j = j0;
  x = i0+.5;
  y = j0+.5;

  // NOTE: Don't need the intermediate storage (x,y). It's just for clarity
  // during development.
  sfore[0][0] = x;
  sfore[0][1] = y;
  dfore[0] = 0.;

  // Take care of first segment separately.
  n = 0;
  ds_x;
  ds_y;
  i = (int)floor( x);
  j = (int)floor( y);

  // Seed pixel.
  tfore[n][0] = i;
  tfore[n][1] = j;

  vx = v[j][i][0];
  vy = v[j][i][1];

  // Tweak vectors that run parallel to cell walls.
  if( fabs(vx) < EPS) { vx+=((vx<0)?(-EPS):(EPS)); vx*=10.;}
  if( fabs(vy) < EPS) { vy+=((vy<0)?(-EPS):(EPS)); vy*=10.;}

  switch( which_quadrant_vector_points_into( vx, vy))
  {
    case 1:
      if( vx > vy) { ds_x = .5; ds_y = ds_x * ( vy / vx);}
      else         { ds_y = .5; ds_x = ds_y * ( vx / vy);}
      ds_x+=EPS; ds_y+=EPS;
      break;
    case 2:
      if( -vx > vy) { ds_x = -.5; ds_y = ds_x * ( vy / vx);}
      else          { ds_y =  .5; ds_x = ds_y * ( vx / vy);}
      ds_x-=EPS; ds_y+=EPS;
      break;
    case 3:
      if( -vx > -vy) { ds_x = -.5; ds_y = ds_x * ( vy / vx);}
      else           { ds_y = -.5; ds_x = ds_y * ( vx / vy);}
      ds_x-=EPS; ds_y-=EPS;
      break;
    case 4:
      if( vx > -vy) { ds_x =  .5; ds_y = ds_x * ( vy / vx);}
      else          { ds_y = -.5; ds_x = ds_y * ( vx / vy);}
      ds_x+=EPS; ds_y-=EPS;
      break;
    default:
      cout << "ERROR: Unhandled case: which_quadrant_vector_points_into("
           << vx << "," << vy << ") = "
           << which_quadrant_vector_points_into( vx, vy) << "." << endl;
      exit(1);
  }

  x += ds_x;
  y += ds_y;
  n++;

  loop = false;

  // Compute the segments in the positive direction.
  while( n <= ns && not_off_domain(x,y) && !loop)
  {
    sfore[n][0] = x;
    sfore[n][1] = y;

    // Adjust dfore for periodicity if periodic==true.
    if( 1 && periodic && fabs( sfore[n][0]-sfore[n-1][0]) > ni/2.)
    {
      if( 0)
      {
        cout << __FILE__ << " " << __LINE__ << " >> "
             << "Crossover from x="
             << sfore[n-1][0] << " to x="
             << sfore[n][0] << "." << endl;
      }
      ds_x = sfore[n][0]+((sfore[n][0]-sfore[n-1][0]<0)?(1.):(-1.))*ni-sfore[n-1][0];
    }
    else
    {
      ds_x = sfore[n][0]-sfore[n-1][0];
    }
    if( 1 && periodic && fabs( sfore[n][1]-sfore[n-1][1]) > nj/2.)
    {
      if( 0)
      {
        cout << __FILE__ << " " << __LINE__ << " >> "
             << "Crossover from y="
             << sfore[n-1][1] << " to y="
             << sfore[n][1] << "." << endl;
      }
      ds_y = sfore[n][1]+((sfore[n][1]-sfore[n-1][1]<0)?(1.):(-1.))*nj-sfore[n-1][1];
    }
    else
    {
      ds_y = sfore[n][1]-sfore[n-1][1];
    }
    dfore[n] = dfore[n-1] + sqrt( (ds_x)*(ds_x) + (ds_y)*(ds_y) );

    n++;

    i = (int)floor( x);
    j = (int)floor( y);

    // Seed pixel.
    tfore[n-1][0] = i;
    tfore[n-1][1] = j;
#if DETECT_LOOPS
    for( k=0; k<n-1; k++)
    {
      if( tfore[k][0]==tfore[n-1][0] && tfore[k][1]==tfore[n-1][1])
      {
        loop=true;
      //cout << "(" << i0 << "," << j0 << ").."
      //     << "(" << tfore[k][0] << "," << tfore[k][1] << "): "
      //     << "loop fore! n = " << n << ", k = " << k << endl;
      }
    }
#endif

    if( !loop)
    {
      if( i<0 || j<0 || i>=ni || j>=nj)
      {
        cout << "(i,j)=(" << i << "," << j << ") is bad!" << endl;
        exit(1);
      }

      vx = v[j][i][0];
      vy = v[j][i][1];

      // Tweak vectors that run parallel to cell walls (to avoid round-off
      // awkwardness).
      if( fabs(vx) < EPS) { vx+=((vx<0)?(-EPS):(EPS)); vx*=10.;}
      if( fabs(vy) < EPS) { vy+=((vy<0)?(-EPS):(EPS)); vy*=10.;}

      switch( which_quadrant_vector_points_into( vx, vy))
      {
        case 1: { // First quadrant.
          if( pointing_at_top_wall( x, y, vx, vy))
          {
            ds_y = ceil(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x+EPS);
            y += ( ds_y+EPS);
          }
          else if( pointing_at_right_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = ceil(x) - x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x+EPS);
            y += ( ds_y+EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        }
        case 2: { // Second quadrant.
          if( pointing_at_top_wall( x, y, vx, vy))
          {
            ds_y = ceil(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x-EPS);
            y += ( ds_y+EPS);
          }
          else if( pointing_at_left_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = floor(x)-x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x-EPS);
            y += ( ds_y+EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        }
        case 3: { // Third quadrant.
          if( pointing_at_bottom_wall( x, y, vx, vy))
          {
            ds_y = floor(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x-EPS);
            y += ( ds_y-EPS);
          }
          else if( pointing_at_left_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = floor(x) - x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x-EPS);
            y += ( ds_y-EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        }
        case 4: { // Fourth quadrant.
          if( pointing_at_bottom_wall( x, y, vx, vy))
          {
            ds_y = floor(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x+EPS);
            y += ( ds_y-EPS);
          }
          else if( pointing_at_right_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = ceil(x) - x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x+EPS);
            y += ( ds_y-EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        }
        default: { // Unhandled case.
          cout << "ERROR: Unhandled case!" << endl;
          exit(1);
        }
      }

    }
  }

  if( not_off_domain( x, y) && !loop)
  {
    i = (int)floor( x);
    j = (int)floor( y);

    // Seed pixel.
    tfore[n-1][0] = i;
    tfore[n-1][1] = j;

    countpos = n;
  }
  else
  {
    countpos = n-1;
    for( ; n<=ns; n++) { dfore[n] = dfore[n-1];}
  }
 }
 else
 {
   countpos = 0;
 }

 if(1)
 {
  //############################################################################
  //
  // B A C K W A R D
  //
  i = i0;
  j = j0;

  x = i0+.5;
  y = j0+.5;
  // NOTE: Don't need the intermediate storage (x,y). It's just for clarity
  // during development.
  sback[0][0] = x;
  sback[0][1] = y;
  dback[0] = 0.;

  // Take care of first segment separately.
  n = 0;
  ds_x;
  ds_y;
  i = (int)floor( x);
  j = (int)floor( y);

  // Seed pixel.
  tback[n][0] = i;
  tback[n][1] = j;

  vx = -v[j][i][0];
  vy = -v[j][i][1];

  // Tweak vectors that run parallel to cell walls.
  if( fabs(vx) < EPS) { vx+=((vx<0)?(-EPS):(EPS)); vx*=10.;}
  if( fabs(vy) < EPS) { vy+=((vy<0)?(-EPS):(EPS)); vy*=10.;}

  switch( which_quadrant_vector_points_into( vx, vy))
  {
    case 1:
      if( vy < vx) { ds_x = .5; ds_y = ds_x * ( vy / vx);}
      else         { ds_y = .5; ds_x = ds_y * ( vx / vy);}
      ds_x+=EPS; ds_y+=EPS;
      break;
    case 2:
      if( vy < -vx) { ds_x = -.5; ds_y = ds_x * ( vy / vx);}
      else          { ds_y =  .5; ds_x = ds_y * ( vx / vy);}
      ds_x-=EPS; ds_y+=EPS;
      break;
    case 3:
      if( -vy < -vx) { ds_x = -.5; ds_y = ds_x * ( vy / vx);}
      else           { ds_y = -.5; ds_x = ds_y * ( vx / vy);}
      ds_x-=EPS; ds_y-=EPS;
      break;
    case 4:
      if( -vy < vx) { ds_x =  .5; ds_y = ds_x * ( vy / vx);}
      else          { ds_y = -.5; ds_x = ds_y * ( vx / vy);}
      ds_x+=EPS; ds_y-=EPS;
      break;
    default:
      cout << "ERROR: Unhandled case: which_quadrant_vector_points_into("
           << vx << "," << vy << ") = "
           << which_quadrant_vector_points_into( vx, vy) << "." << endl;
      exit(1);
  }

  x += ds_x;
  y += ds_y;
  n++;

  loop = false;

  // Compute the segments in the positive direction.
  while( n <= ns && not_off_domain(x,y) && !loop)
  {
    sback[n][0] = x;
    sback[n][1] = y;

    // Adjust dback for periodicity if periodic==true.
    if( 1 && periodic && fabs( sback[n][0]-sback[n-1][0]) > ni/2.)
    {
      ds_x = sback[n][0]+((sback[n][0]-sback[n-1][0]<0)?(1.):(-1.))*ni-sback[n-1][0];
    }
    else
    {
      ds_x = sback[n][0]-sback[n-1][0];
    }
    if( 1 && periodic && fabs( sback[n][1]-sback[n-1][1]) > nj/2.)
    {
      ds_y = sback[n][1]+((sback[n][1]-sback[n-1][1]<0)?(1.):(-1.))*nj-sback[n-1][1];
    }
    else
    {
      ds_y = sback[n][1]-sback[n-1][1];
    }
    dback[n] = dback[n-1] + sqrt( (ds_x)*(ds_x) + (ds_y)*(ds_y) );

    n++;

    i = (int)floor( x);
    j = (int)floor( y);

    // Seed pixel.
    tback[n-1][0] = i;
    tback[n-1][1] = j;
#if DETECT_LOOPS
    for( k=0; k<n-1; k++)
    {
      if( tback[k][0]==tback[n-1][0] && tback[k][1]==tback[n-1][1])
      {
        loop=true;
      //cout << "(" << i0 << "," << j0 << ").."
      //     << "(" << tback[k][0] << "," << tback[k][1] << "): "
      //     << "loop back! n = " << n << ", k = " << k << endl;
      }
    }
#endif

    if( !loop)
    {
      vx = -v[j][i][0];
      vy = -v[j][i][1];

      // Tweak vectors that run parallel to cell walls.
      if( fabs(vx) < EPS) { vx+=((vx<0)?(-EPS):(EPS)); vx*=10.;}
      if( fabs(vy) < EPS) { vy+=((vy<0)?(-EPS):(EPS)); vy*=10.;}

      switch( which_quadrant_vector_points_into( vx, vy))
      {
        case 1: // First quadrant.
          if( pointing_at_top_wall( x, y, vx, vy))
          {
            ds_y = ceil(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x+EPS);
            y += ( ds_y+EPS);
          }
          else if( pointing_at_right_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = ceil(x) - x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x+EPS);
            y += ( ds_y+EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        case 2: // Second quadrant.
          if( pointing_at_top_wall( x, y, vx, vy))
          {
            ds_y = ceil(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x-EPS);
            y += ( ds_y+EPS);
          }
          else if( pointing_at_left_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = floor(x)-x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x-EPS);
            y += ( ds_y+EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        case 3: // Third quadrant.
          if( pointing_at_bottom_wall( x, y, vx, vy))
          {
            ds_y = floor(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x-EPS);
            y += ( ds_y-EPS);
          }
          else if( pointing_at_left_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = floor(x) - x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x-EPS);
            y += ( ds_y-EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        case 4: // Fourth quadrant.
          if( pointing_at_bottom_wall( x, y, vx, vy))
          {
            ds_y = floor(y) - y;
            ds_x = ds_y * ( vx / vy);
            x += ( ds_x+EPS);
            y += ( ds_y-EPS);
          }
          else if( pointing_at_right_wall( x, y, vx, vy)) // Can omit the 'if' here.
          {
            ds_x = ceil(x) - x;
            ds_y = ds_x * ( vy / vx);
            x += ( ds_x+EPS);
            y += ( ds_y-EPS);
          }
          else
          {
            cout << "ERROR: This shouldn't happen!" << endl;
            exit(1);
          }
          break;
        default:
          cout << "ERROR: Unhandled case!" << endl;
          exit(1);

      } // switch( which_quadrant_vector_points_into( vx, vy))
    } // if( !loop)
  } // while( n <= ns && not_off_domain(x,y) && !loop)

  if( not_off_domain( x, y) && !loop)
  {
    i = (int)floor( x);
    j = (int)floor( y);

    // Seed pixel.
    tback[n-1][0] = i;
    tback[n-1][1] = j;

    countneg = n;
  }
  else
  {
    countneg = n-1;
    for( ; n<=ns; n++) { dback[n] = dback[n-1];}
  }
 }
 else
 {
   countneg = 0;
 }

} // void lic_class::compute_streamline(

lic_class::value_type lic_class::compute_h_fore_approx(
  const int n0,
  const int countpos)
{
  if( n0 < countpos)
  {
    if( tfore[n0-1][0] == tfore[n0][0])
    {
      if( tfore[n0-1][0] == tfore[n0+1][0])
      {
        // Straight through
        return 1.;
      }
      else
      {
        // Corner
        return 1./sqrt(2.);
      }
    }
    else
    {
      if( tfore[n0-1][1] == tfore[n0+1][1])
      {
        // Straight through
        return 1.;
      }
      else
      {
        // Corner
        return 1./sqrt(2.);
      }
    }
  }
  else // Last pixel, so assume straight through.
  {
    return 1.;
  }
}

lic_class::value_type lic_class::compute_h_back_approx(
  const int n0,
  const int countneg)
{
  if( n0 < countneg)
  {
    if( tback[n0-1][0] == tback[n0][0])
    {
      if( tback[n0-1][0] == tback[n0+1][0])
      {
        // Straight through
        return 1.;
      }
      else
      {
        // Corner
        return 1./sqrt(2.);
      }
    }
    else
    {
      if( tback[n0-1][1] == tback[n0+1][1])
      {
        // Straight through
        return 1.;
      }
      else
      {
        // Corner
        return 1./sqrt(2.);
      }
    }
  }
  else // Last pixel, so assume straight through.
  {
    return 1.;
  }

}

lic_class::value_type lic_class::compute_h_fore( const int n)
{
  //return 1.;
#if 0
  return sqrt( (sfore[n+1][0]-sfore[n][0])*(sfore[n+1][0]-sfore[n][0])
             + (sfore[n+1][1]-sfore[n][1])*(sfore[n+1][1]-sfore[n][1]) );
#else
#if 0
  return
    sqrt( (sfore[n+1][0]-sfore[n][0])*(sfore[n+1][0]-sfore[n][0])
        + (sfore[n+1][1]-sfore[n][1])*(sfore[n+1][1]-sfore[n][1]) )
    *
    .5*( ((dfore[n+1]<20.)?(20.-dfore[n+1]):(0.))
       + ((dfore[n  ]<20.)?(20.-dfore[n  ]):(0.)))
    //.5*((dfore[ns]-dfore[n+1])+(dfore[ns]-dfore[n]))
    //(sfore[n+1][0]-sfore[n][0]) (sfore[n+1][1]-sfore[n][1])
  ;
#else
//cout << ".5*( (dfore[ns]-dfore[n+1]) + (dfore[ns]-dfore[n])) = "
//     << ".5*( (" << dfore[ns] << "-" << dfore[n+1] << ") + ("
//                 << dfore[ns] << "-" << dfore[n] << ")) = "
//     << .5*( (dfore[ns]-dfore[n+1]) + (dfore[ns]-dfore[n])) << endl;

  // Adjust for periodicity if periodic==true.
  real ds_x, ds_y;
  if( 1 && periodic && fabs( sfore[n+1][0]-sfore[n][0]) > ni/2.)
  {
    ds_x = sfore[n+1][0]+((sfore[n+1][0]-sfore[n][0]<0)?(1.):(-1.))*ni-sfore[n][0];
  }
  else
  {
    ds_x = sfore[n+1][0]-sfore[n][0];
  }
  if( 1 && periodic && fabs( sfore[n+1][1]-sfore[n][1]) > nj/2.)
  {
    ds_y = sfore[n+1][1]+((sfore[n+1][1]-sfore[n][1]<0)?(1.):(-1.))*nj-sfore[n][1];
  }
  else
  {
    ds_y = sfore[n+1][1]-sfore[n][1];
  }
  return
    sqrt( ds_x*ds_x + ds_y*ds_y )
   *
    .5*( (dfore[ns]-dfore[n+1]) + (dfore[ns]-dfore[n]));
#endif
#endif
}

lic_class::value_type lic_class::compute_h_back( const int n)
{
  //return 1.;
#if 0
  return sqrt( (sback[n+1][0]-sback[n][0])*(sback[n+1][0]-sback[n][0])
             + (sback[n+1][1]-sback[n][1])*(sback[n+1][1]-sback[n][1]) );
#else
#if 0
  return
    sqrt( (sback[n+1][0]-sback[n][0])*(sback[n+1][0]-sback[n][0])
        + (sback[n+1][1]-sback[n][1])*(sback[n+1][1]-sback[n][1]) )
    *
    .5*( ((dback[n+1]<20.)?(20.-dback[n+1]):(0.))
       + ((dback[n  ]<20.)?(20.-dback[n  ]):(0.)))
    //.5*((dback[ns]-dback[n+1])+(dback[ns]-dback[n]))
    //(sback[n+1][0]-sback[n][0]) (sback[n+1][1]-sback[n][1])
  ;
#else

  // Adjust for periodicity if periodic==true.
  real ds_x, ds_y;
  if( 1 && periodic && fabs( sback[n+1][0]-sback[n][0]) > ni/2.)
  {
    ds_x = sback[n+1][0]+((sback[n+1][0]-sback[n][0]<0)?(1.):(-1.))*ni-sback[n][0];
  }
  else
  {
    ds_x = sback[n+1][0]-sback[n][0];
  }
  if( 1 && periodic && fabs( sback[n+1][1]-sback[n][1]) > nj/2.)
  {
    ds_y = sback[n+1][1]+((sback[n+1][1]-sback[n][1]<0)?(1.):(-1.))*nj-sback[n][1];
  }
  else
  {
    ds_y = sback[n+1][1]-sback[n][1];
  }

  return
    sqrt( ds_x*ds_x + ds_y*ds_y )
   *
    .5*( (dback[ns]-dback[n+1]) + (dback[ns]-dback[n]));
#endif
#endif
}

bool lic_class::are_valid_indices(
       const int i, const int j) const
{
  return ( i>=0 && i<ni && j>=0 && j<nj);
}

bool lic_class::pointing_at_top_wall(
  const real x,
  const real y,
  const real vx,
  const real vy ) const
{
  switch( which_quadrant_vector_points_into( vx, vy))
  {
    case 1:
      if( (ceil(x)-x)*(vy/vx) < (ceil(y)-y)) { return false;}
                                        else { return true;}
      break;
    case 2:
      if( (floor(x)-x)*(vy/vx) < (ceil(y)-y)) { return false;}
                                         else { return true;}
      break;
    case 3: return false; break;
    case 4: return false; break;
    default:
      cout << "ERROR: Unhandled case!" << endl;
      exit(1);
  }
}

bool lic_class::pointing_at_bottom_wall(
  const real x,
  const real y,
  const real vx,
  const real vy ) const
{
  switch( which_quadrant_vector_points_into( vx, vy))
  {
    case 1: return false; break;
    case 2: return false; break;
    case 3:
      if( (x-floor(x))*(vy/vx) < (y-floor(y))) { return false;}
                                          else { return true;}
      break;
    case 4:
      if( (ceil(x)-x)*(-vy/vx) < (y-floor(y))) { return false;}
                                          else { return true;}
      break;
    default:
      cout << "ERROR: Unhandled case!" << endl;
      exit(1);
  }
}

bool lic_class::pointing_at_right_wall(
  const real x,
  const real y,
  const real vx,
  const real vy ) const
{
  switch( which_quadrant_vector_points_into( vx, vy))
  {
    case 1:
      if( (ceil(x)-x)*(vy/vx) < (ceil(y)-y)) { return true;}
                                           else { return false;}
      break;
    case 2: return false; break;
    case 3: return false; break;
    case 4:
      if( (ceil(x)-x)*(-vy/vx) < (y-floor(y))) { return true;}
                                            else { return false;}
      break;
    default:
      cout << "ERROR: Unhandled case!" << endl;
      exit(1);
  }
}

bool lic_class::pointing_at_left_wall(
  const real x,
  const real y,
  const real vx,
  const real vy ) const
{
  switch( which_quadrant_vector_points_into( vx, vy))
  {
    case 1: return false; break;
    case 2:
      if( (floor(x)-x)*(vy/vx) < (ceil(y)-y)) { return true;}
                                               else { return false;}
      break;
    case 3:
      if( (floor(x)-x)*(-vy/vx) < (y-floor(y))) { return true;}
                                                else { return false;}
      break;
    case 4: return false; break;
    default:
      cout << "ERROR: Unhandled case!" << endl;
      exit(1);
  }
}


bool lic_class::not_off_domain( real &x, real &y) const
{
  if( !periodic)
  {
    return ( x > 0 && x < ni && y > 0 && y < nj);
  }
  else
  {
    x = (x<0)?(ni+x):(x);
    x = (x>ni)?(x-ni):(x);
    y = (y<0)?(nj+y):(y);
    y = (y>nj)?(y-nj):(y);
    return true;
  }
}

int lic_class::which_quadrant_vector_points_into(
      const real vx,
      const real vy ) const
{
  if( vx > 0) { if( vy > 0) { return 1;} else { return 4;}}
  else        { if( vy > 0) { return 2;} else { return 3;}}
}

lic_class::~lic_class()
{
  int i, j;

  // Free v.
  for( j=0; j<nj; j++)
  {
    for( i=0; i<ni; i++)
    {
      delete [] v[j][i];
    }
    delete [] v[j];
  }
  delete [] v;

  // Free r.
  for( j=0; j<nj; j++)
  {
    delete [] r[j];
  }
  delete [] r;

  // Free rprime, gprime, bprime.
  for( j=0; j<nj; j++)
  {
    delete [] rprime[j];
    delete [] gprime[j];
    delete [] bprime[j];
  }
  delete [] rprime;
  delete [] gprime;
  delete [] bprime;

  // Free tfore and tback.
  for( i=0; i<ns; i++)
  {
    delete [] tfore[i];
    delete [] tback[i];
  }
  delete [] tfore;
  delete [] tback;

  // Free sfore and sback.
  for( i=0; i<ns; i++)
  {
    delete [] sfore[i];
    delete [] sback[i];
  }
  delete [] sfore;
  delete [] sback;

  // Free dfore and dback.
  delete [] dfore;
  delete [] dback;

}

ostream &operator<<( ostream &o, const lic_class &arg_lic)
{
  arg_lic.display( o);
  return o;
}

} //namespace dthorne0_lic_class_ns

// vim: foldmethod=syntax foldlevel=1

//##############################################################################
//
// lic.cc
//

#include "lic_class.h"
using namespace dthorne0_lic_class_ns;

#include <iostream>
using namespace std;

int main( int argc, char **argv)
{
  cout << endl << "main(): Hello!" << endl << endl;

  if( argc != 1)
  {
    switch( argc)
    {
      case 2:
      {
        if( isdigit(argv[1][0]))
        {
          cout << "argv[1] = " << argv[1] << endl;
          cout << "atoi(argv[1]) = " << atoi(argv[1]) << endl;
          lic_class lic( (lic_class::texture_type)atoi(argv[1]));
          cout << lic << endl;
        }
        else
        {
          lic_class lic( 3, argv[1]);
          cout << lic << endl;
        }
        break;
      }
      case 3:
      {
        if( isdigit(argv[1][0]))
        {
          lic_class lic( 3, argv[1], atoi(argv[2]));
          cout << lic << endl;
        }
        else if( isdigit(argv[2][0]))
        {
          lic_class lic( atoi(argv[2]), argv[1]);
          cout << lic << endl;
        }
        else
        {
          lic_class lic( 3, argv[1]);
          cout << lic << endl;
        }
        break;
      }
      case 4:
      {
        lic_class lic( atoi(argv[1]), argv[2], atoi(argv[3]));
        cout << lic << endl;
        break;
      }
      default:
      {
        cout << "usage: " << argv[0] << " [filename]" << endl;
        exit(1);
      }
    }
  }
  else
  {
    lic_class lic;
    cout << lic << endl;
  }

  cout << endl << "main(): Terminating normally!" << endl << endl;

  return 0;
}

// vim: foldmethod=syntax foldlevel=1

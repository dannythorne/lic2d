//##############################################################################
//
// bing.h
//

#ifndef BING_H
#define BING_H

#define BING( str)\
  std::cout\
    << __FILE__\
    << " "\
    << __LINE__\
    << " >> "\
    << str\
    << "\n";

#if DEBUG > 0
#define BING1( str)\
  std::cout\
    << __FILE__\
    << " "\
    << __LINE__\
    << " >> "\
    << str\
    << "\n";
#else
#define BING1( str)
#endif

#if DEBUG > 1
#define BING2( str)\
  std::cout\
    << "  " \
    << __FILE__\
    << " "\
    << __LINE__\
    << " >> "\
    << str\
    << "\n";
#else
#define BING2( str)
#endif

#endif

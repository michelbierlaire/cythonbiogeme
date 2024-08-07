//-*-c++-*------------------------------------------------------------
//
// File name : bioDebug.h
// @date   Wed Apr 11 13:19:50 2018
// @author Michel Bierlaire
// @version Revision 1.0
//
//--------------------------------------------------------------------

    #ifndef bioDebug_h
    #define bioDebug_h

    #include <iostream>

    #ifdef DEBUG
    #define DEBUG_MESSAGE(message) do { std::cout << __FILE__ << ":" << __LINE__ << " " << message << std::endl; } while (0)
    #else
    #define DEBUG_MESSAGE(message) do {} while (0)
    #endif

    #endif

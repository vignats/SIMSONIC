#ifndef _COUTWINDOWS_HPP_
#define _COUTWINDOWS_HPP_

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
#include <iostream>
#include "mex.h"

class mystream : public std::streambuf
{
protected:
    virtual std::streamsize xsputn(const char *s, std::streamsize n)
    {
        mexPrintf("%.*s", n, s);
        return n;
    }
    virtual int overflow(int c = EOF)
    {
        if (c != EOF)
        {
            mexPrintf("%.1s", &c);
        }
        return 1;
    }
};
class scoped_redirect_cout
{
public:
    scoped_redirect_cout()
    {
        old_buf = std::cout.rdbuf();
        std::cout.rdbuf(&mout);
    }
    ~scoped_redirect_cout() { std::cout.rdbuf(old_buf); }

private:
    mystream mout;
    std::streambuf *old_buf;
};
static scoped_redirect_cout mycout_redirect;

#endif
#endif

//Tokenizer.h

#ifndef __TOKENIZER_H__
#define __TOKENIZER_H__
#include <iostream>
#include <functional>
#include <string>
#include <vector>
#include <algorithm>
#include <locale>

//For the case the default is a space.
//This is the default predicate for the Tokenize() function.
class CIsSpace : public std::unary_function<char, bool>
{
public:
    bool operator()(char c) const;
};

inline bool CIsSpace::operator()(char c) const
{
    //isspace<char> returns true if c is a white-space character (0x09-0x0D or 0x20)
    //return isspace<char>(c);
    return isspace(c);
}

//For the case the separator is a comma
class CIsComma : public std::unary_function<char, bool>
{
public:
    bool operator()(char c) const;
};

inline bool CIsComma::operator()(char c) const
{
    return (',' == c);
}

//For the case the separator is a colon
class CIsColon : public std::unary_function<char, bool>
{
public:
    bool operator()(char c) const;
};

inline bool CIsColon::operator()(char c) const
{
    return (':' == c);
}

//For the case the separator is a colon
class CIsSemiColon : public std::unary_function<char, bool>
{
public:
    bool operator()(char c) const;
};

inline bool CIsSemiColon::operator()(char c) const
{
    return (';' == c);
}



//For the case the separator is a character from a set of characters given in a string
class CIsFromString : public std::unary_function<char, bool>
{
public:
    CIsFromString(std::string const &rostr) : m_ostr(rostr) {}
    bool operator()(char c) const;
    
private:
    std::string m_ostr;
};

inline bool CIsFromString::operator()(char c) const
{
    std::size_t iFind = m_ostr.find(c);
    if(iFind != std::string::npos)
        return true;
    else
        return false;
}

//String Tokenizer
template <class Pred = CIsSpace>
class CTokenizer
{
    public:
    //The predicate should evaluate to true when applied to a separator.
    static void Tokenize(std::vector<std::string>& roResult, std::string const& rostr, Pred const& roPred=Pred() );
};

//The predicate should evaluate to true when applied to a separator.
template <class Pred>
inline void CTokenizer<Pred>::Tokenize(std::vector<std::string>& roResult, std::string const& rostr, Pred const& roPred)
{
    //First clear the results vector
    roResult.clear();
    std::string::const_iterator it = rostr.begin();
    std::string::const_iterator itTokenEnd = rostr.begin();
    while(it != rostr.end() ) {
        //Eat seperators
        while(roPred(*it))
            it++;
        //Find next token
        itTokenEnd = find_if(it, rostr.end(), roPred);
        //Append token to result
        if(it < itTokenEnd)
            roResult.push_back(std::string(it, itTokenEnd));
        it = itTokenEnd;
    }
}
#endif //__TOKENIZER_H__





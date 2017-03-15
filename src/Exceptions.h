#ifndef EXCEPTIONS_H_
#define EXCEPTIONS_H_

#include <stdexcept>
#include <string>


struct IntegrationException: public std::runtime_error
{
    IntegrationException(std::string const& message)
        : std::runtime_error(message + " Was thrown")
    {}
};


#endif /* EXCEPTIONS_H_ */

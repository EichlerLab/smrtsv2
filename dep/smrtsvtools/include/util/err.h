/*
 * Header for error reporting utilities
 */

#ifndef INCLUDE_UTIL_ERR_H_
#define INCLUDE_UTIL_ERR_H_

/**
 * Print a formatted error message to STDERR.
 *
 * @param message Message with formatting characters.
 * @param ... Format arguments.
 */
void err(const char *message, ...);

/**
 * Print a formatted error message to STDERR.
 *
 * @param message Message with formatting characters.
 */
void err(const std::exception &ex);

/**
 * Print a formatted warning message to STDERR.
 *
 * @param message Message with formatting characters.
 * @param ... Format arguments.
 */
void warn(const char *message, ...);

/**
 * Print a formatted warning message to STDERR.
 *
 * @param message Message with formatting characters.
 */
void warn(const std::exception &ex);



#endif /* INCLUDE_UTIL_ERR_H_ */

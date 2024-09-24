/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#pragma once

#ifdef GOTOOLS_LOG
#include <spdlog/spdlog.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/details/registry.h>
#include "spdlog/fmt/fmt.h"
//#include <fmt/core.h>
//#include <fmt/format.h>
#include <iostream> // Added this line
#include <memory> // Include for std::shared_ptr
#include <fstream>
#else
#include <iostream>
#endif

namespace Go {

#ifdef GOTOOLS_LOG

class Logger {
public:
    static void init(const std::string& logfile_name = "logfile.txt") {
        static bool initialized = false; // Track initialization status
        if (!initialized) {
            try {
                std::cout << "Setting logfile: " << logfile_name << std::endl;
                // Open the file in truncate mode to clear its contents
                file_logger = spdlog::rotating_logger_mt("file_logger", logfile_name, 1048576 * 5, 3); // 5 MB size, 3 files
                std::ofstream(logfile_name, std::ios::trunc).close(); // Truncate the file
                spdlog::set_default_logger(file_logger);
                spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
                spdlog::set_level(spdlog::level::warn); // Default log level
                initialized = true; // Mark as initialized
            } catch (const spdlog::spdlog_ex& ex) {
                std::cerr << "Log initialization failed: " << ex.what() << std::endl;
            }
        }// else {
        //     std::cout << "Logger has already been initialized." << std::endl; // Optional: log a message if already initialized
        //}
    }

    // Enum log levels (0-6): trace, debug, info, warn, err, critical, off.
    static void setLogLevel(size_t level) {
        spdlog::level::level_enum log_level = static_cast<spdlog::level::level_enum>(level); // Convert int to level_enum
        spdlog::set_level(log_level);
        std::cout << "Log level set to: " << level << std::endl; // Explicitly cast level to int
    }

    // Logger pointer
    static std::shared_ptr<spdlog::logger> file_logger; // Declare the logger

private:
    Logger() = default; // Prevent instantiation

    // Static instance to ensure automatic initialization
    struct LoggerInitializer {
        LoggerInitializer() {
            init(); // Call init() automatically
        }
    };

    static LoggerInitializer loggerInitializer; // Create a static instance of LoggerInitializer
};

// Move macro definitions below the Logger class definition

#define LOG_TRACE(...)   \
    do {                 \
        Go::Logger::init();  \
        spdlog::trace(__VA_ARGS__); \
    } while (0)

#define LOG_DEBUG(...)   \
    do {                 \
        Go::Logger::init();  \
        spdlog::debug(__VA_ARGS__); \
    } while (0)

#define LOG_INFO(...)    \
    do {                 \
        Go::Logger::init();  \
        spdlog::info(__VA_ARGS__); \
    } while (0)

#define LOG_WARN(...)    \
    do {                 \
        Go::Logger::init();  \
        spdlog::warn(__VA_ARGS__); \
    } while (0)

#define LOG_ERROR(...)   \
    do {                 \
        Go::Logger::init();  \
        spdlog::error(__VA_ARGS__); \
    } while (0)

#define LOG_CRITICAL(...) \
    do {                  \
        Go::Logger::init();   \
        spdlog::critical(__VA_ARGS__); \
    } while (0)

#else

class Logger {
public:
    static void init(const std::string& logfile_name = "logfile.txt")
    {
        std::cout << "Logging not enabled. Sending all log messages to cerr." << std::endl;
    }
    
    // Enum log levels (0-6): trace, debug, info, warn, err, critical, off.
    static void setLogLevel(size_t level) {
        log_level = level; // Set the log level
        std::cout << "Log level set to: " << log_level << std::endl; // Use the original int value
    }

    static size_t log_level;

private:
    Logger() = default; // Prevent instantiation

};

#define LOG_TRACE(...)   if (Go::Logger::log_level <= 0) std::cerr << "[TRACE] " << __VA_ARGS__ << std::endl
#define LOG_DEBUG(...)   if (Go::Logger::log_level <= 1) std::cerr << "[DEBUG] " << __VA_ARGS__ << std::endl
#define LOG_INFO(...)    if (Go::Logger::log_level <= 2) std::cerr << "[INFO] " << __VA_ARGS__ << std::endl
#define LOG_WARN(...)    if (Go::Logger::log_level <= 3) std::cerr << "[WARN] " << __VA_ARGS__ << std::endl
#define LOG_ERROR(...)   if (Go::Logger::log_level <= 4) std::cerr << "[ERROR] " << __VA_ARGS__ << std::endl
#define LOG_CRITICAL(...) if (Go::Logger::log_level <= 5) std::cerr << "[CRITICAL] " << __VA_ARGS__ << std::endl

#endif

} // namespace Go
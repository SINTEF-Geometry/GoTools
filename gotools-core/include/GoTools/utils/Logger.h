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

class Logger {
public:
    static void init(const std::string& logfile_name = "logfile.txt") {
        static bool initialized = false; // Track initialization status
        if (!initialized) {
            try {
                // Open the file in truncate mode to clear its contents
                file_logger = spdlog::rotating_logger_mt("file_logger", logfile_name, 1048576 * 5, 3); // 5 MB size, 3 files
                std::ofstream(logfile_name, std::ios::trunc).close(); // Truncate the file
                spdlog::set_default_logger(file_logger);
                spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
                spdlog::set_level(spdlog::level::trace); // Default log level
                initialized = true; // Mark as initialized
            } catch (const spdlog::spdlog_ex& ex) {
                std::cerr << "Log initialization failed: " << ex.what() << std::endl;
            }
        }// else {
        //     std::cout << "Logger has already been initialized." << std::endl; // Optional: log a message if already initialized
        //}
    }

    static void setLogLevel(spdlog::level::level_enum level) {
        spdlog::set_level(level);
        std::cout << "Log level set to: " << static_cast<int>(level) << std::endl; // Explicitly cast level to int
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
#ifdef GOTOOLS_LOG

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

inline void init(const std::string& logfile_name = "logfile.txt") {} // Empty init function when logging is disabled
inline void setLogLevel(spdlog::level::level_enum) {} // Empty setLogLevel function when logging is disabled

#define LOG_TRACE(...)   std::cout << "[TRACE] " << __VA_ARGS__ << std::endl
#define LOG_DEBUG(...)   std::cout << "[DEBUG] " << __VA_ARGS__ << std::endl
#define LOG_INFO(...)    std::cout << "[INFO] " << __VA_ARGS__ << std::endl
#define LOG_WARN(...)    std::cerr << "[WARN] " << __VA_ARGS__ << std::endl
#define LOG_ERROR(...)   std::cerr << "[ERROR] " << __VA_ARGS__ << std::endl
#define LOG_CRITICAL(...) std::cerr << "[CRITICAL] " << __VA_ARGS__ << std::endl

#endif

} // namespace Go
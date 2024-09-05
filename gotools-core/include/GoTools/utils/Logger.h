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

#else
#include <iostream>
#endif

namespace Go::Logger {

#ifdef GOTOOLS_LOG

// Declare a shared pointer to hold the logger within the Go::Logger namespace
extern std::shared_ptr<spdlog::logger> file_logger; // Change to extern

inline void init(const std::string& logfile_name = "logfile.txt") {
    static bool initialized = false;
    if (!initialized) {
        try {
            // Create a rotating file sink
            auto file_logger = spdlog::rotating_logger_mt("file_logger", logfile_name, 1048576 * 5, 3); // 5 MB size, 3 files
            spdlog::set_default_logger(file_logger);
            spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
            spdlog::set_level(spdlog::level::trace); // Default log level
            initialized = true;
        } catch (const spdlog::spdlog_ex& ex) {
            std::cerr << "Log initialization failed: " << ex.what() << std::endl;
        }
    }
}

inline void setLogFileName(const std::string& logfile_name) {
    if (Go::Logger::file_logger != nullptr) {
        // Store the current log level
        auto current_level = Go::Logger::file_logger->level(); // Get the current log level from the existing logger
        file_logger->flush(); // Ensure all logs are flushed before changing
        spdlog::drop("file_logger"); // Remove the existing logger
        file_logger = spdlog::rotating_logger_mt("file_logger", logfile_name, 1048576 * 5, 3);
        file_logger->set_level(current_level); // Restore the previous log level
    } else {
        // Create a new logger with the new logfile name
        file_logger = spdlog::rotating_logger_mt("file_logger", logfile_name, 1048576 * 5, 3);
    }
}

inline void setLogLevel(spdlog::level::level_enum level) {
    spdlog::set_level(level);
    std::cout << "Log level set to: " << static_cast<int>(level) << std::endl; // Explicitly cast level to int
}

struct LoggerInitializer {
    LoggerInitializer() {
        init(); // Call init() automatically
    }
};

// Create a static instance of LoggerInitializer
static LoggerInitializer loggerInitializer;

#define LOG_TRACE spdlog::trace
#define LOG_DEBUG spdlog::debug
#define LOG_INFO spdlog::info
#define LOG_WARN spdlog::warn
#define LOG_ERROR spdlog::error
#define LOG_CRITICAL spdlog::critical

#else

inline void init(const std::string& logfile_name = "logfile.txt") {} // Empty init function when logging is disabled
inline void setLogFileName(const std::string&) {} // Empty setLogFileName function when logging is disabled
inline void setLogLevel(spdlog::level::level_enum) {} // Empty setLogLevel function when logging is disabled

#define LOG_TRACE(...)   std::cout << "[TRACE] " << __VA_ARGS__ << std::endl
#define LOG_DEBUG(...)   std::cout << "[DEBUG] " << __VA_ARGS__ << std::endl
#define LOG_INFO(...)    std::cout << "[INFO] " << __VA_ARGS__ << std::endl
#define LOG_WARN(...)    std::cerr << "[WARN] " << __VA_ARGS__ << std::endl
#define LOG_ERROR(...)   std::cerr << "[ERROR] " << __VA_ARGS__ << std::endl
#define LOG_CRITICAL(...) std::cerr << "[CRITICAL] " << __VA_ARGS__ << std::endl

#endif

} // namespace Go::Logger
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
#include "spdlog/sinks/basic_file_sink.h"
#else
#include <iostream>
#endif

namespace Go::Logger {
   
#ifdef GOTOOLS_LOG

inline void init() {
        static bool initialized = false;
        if (!initialized) {
            try {
                auto file_logger = spdlog::basic_logger_mt("file_logger", "logfile.txt");
                spdlog::set_default_logger(file_logger);
                spdlog::set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] %v");
                spdlog::set_level(spdlog::level::trace);
                initialized = true;
                std::cout << "Log initialization succeeded" << std::endl;
            } catch (const spdlog::spdlog_ex& ex) {
                std::cerr << "Log initialization failed: " << ex.what() << std::endl;
            }
        }
    }
 
    #define LOG_TRACE(...) (Go::Logger::init(), SPDLOG_TRACE(__VA_ARGS__))
    #define LOG_DEBUG(...) (Go::Logger::init(), SPDLOG_DEBUG(__VA_ARGS__))
    #define LOG_INFO(...) (Go::Logger::init(), SPDLOG_INFO(__VA_ARGS__))
    #define LOG_WARN(...) (Go::Logger::init(), SPDLOG_WARN(__VA_ARGS__))
    #define LOG_ERROR(...) (Go::Logger::init(), SPDLOG_ERROR(__VA_ARGS__))
    #define LOG_CRITICAL(...) (Go::Logger::init(), SPDLOG_CRITICAL(__VA_ARGS__))
#else
    #define LOG_TRACE(...)   std::cout << "[TRACE] " << __VA_ARGS__ << std::endl
    #define LOG_DEBUG(...)   std::cout << "[DEBUG] " << __VA_ARGS__ << std::endl
    #define LOG_INFO(...)    std::cout << "[INFO] " << __VA_ARGS__ << std::endl
    #define LOG_WARN(...)    std::cerr << "[WARN] " << __VA_ARGS__ << std::endl
    #define LOG_ERROR(...)   std::cerr << "[ERROR] " << __VA_ARGS__ << std::endl
    #define LOG_CRITICAL(...) std::cerr << "[CRITICAL] " << __VA_ARGS__ << std::endl
#endif

    void init();
}
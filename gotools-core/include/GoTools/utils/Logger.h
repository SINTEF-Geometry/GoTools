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
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/support/date_time.hpp>
#else
#include <iostream>
#endif

namespace Go::Logger {
   
#ifdef GOTOOLS_LOG

namespace logging = boost::log;
namespace src = boost::log::sources;
namespace expr = boost::log::expressions;
namespace keywords = boost::log::keywords;

inline void init() {
    static bool initialized = false;
    if (!initialized) {
        try {
            logging::add_file_log(
                keywords::file_name = "logfile.txt",
                keywords::format = (
                    expr::stream
                        << expr::format_date_time< boost::posix_time::ptime >("TimeStamp", "[%Y-%m-%d %H:%M:%S.%f]")
                        << " [" << logging::trivial::severity << "] "
                        << expr::smessage
                )
            );

            logging::core::get()->set_filter(
                logging::trivial::severity >= logging::trivial::trace
            );

            logging::add_common_attributes();

            initialized = true;
            std::cout << "Log initialization succeeded" << std::endl;
        } catch (const std::exception& ex) {
            std::cerr << "Log initialization failed: " << ex.what() << std::endl;
        }
    }
}

#define LOG_TRACE BOOST_LOG_TRIVIAL(trace)
#define LOG_DEBUG BOOST_LOG_TRIVIAL(debug)
#define LOG_INFO BOOST_LOG_TRIVIAL(info)
#define LOG_WARN BOOST_LOG_TRIVIAL(warning)
#define LOG_ERROR BOOST_LOG_TRIVIAL(error)
#define LOG_CRITICAL BOOST_LOG_TRIVIAL(fatal)

#else

inline void init() {} // Empty init function when logging is disabled

#define LOG_TRACE(...)   std::cout << "[TRACE] " << __VA_ARGS__ << std::endl
#define LOG_DEBUG(...)   std::cout << "[DEBUG] " << __VA_ARGS__ << std::endl
#define LOG_INFO(...)    std::cout << "[INFO] " << __VA_ARGS__ << std::endl
#define LOG_WARN(...)    std::cerr << "[WARN] " << __VA_ARGS__ << std::endl
#define LOG_ERROR(...)   std::cerr << "[ERROR] " << __VA_ARGS__ << std::endl
#define LOG_CRITICAL(...) std::cerr << "[CRITICAL] " << __VA_ARGS__ << std::endl
#endif

} // namespace Go::Logger
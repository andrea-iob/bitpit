/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2016 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitbit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
#ifndef __BITPIT_BINARY_UTILS_HPP__
#define __BITPIT_BINARY_UTILS_HPP__

#include <array>
#include <iostream>
#include <limits>
#include <fstream>

namespace bitpit {

/*!
    \brief The namespace 'binary' contains routines for handling binary
    archives.
*/
namespace binary {

    template<typename T>
    void write(std::ostream &stream, const T &value);

    template<typename T>
    void write(std::ostream &stream, const T &value, size_t size);

    template<typename T>
    void write(std::ostream &stream, const T *value);

    template<typename T>
    void write(std::ostream &stream, const T *value, size_t size);

    void write(std::ostream &stream, const std::string &string);

    template<typename T>
    void read(std::istream &stream, T &value);

    template<typename T>
    void read(std::istream &stream, T &value, size_t size);

    template<typename T>
    void read(std::istream &stream, T *value);

    template<typename T>
    void read(std::istream &stream, T *value, size_t size);

    void read(std::istream &stream, std::string &string);

}

}

// Include templates' implementation
#include "binary_utils.tpp"

#endif

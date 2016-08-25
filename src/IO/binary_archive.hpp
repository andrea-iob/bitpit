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
#ifndef __BITPIT_BINARY_ARCHIVE_HPP__
#define __BITPIT_BINARY_ARCHIVE_HPP__

#include <array>
#include <iostream>
#include <limits>
#include <fstream>

namespace bitpit {

class BinaryArchive : protected std::fstream
{

public:
    static const int HEADER_SIZE = 1024;

    static const int VERSION_UNKNOW = - std::numeric_limits<int>::max();

    using std::fstream::close;

    BinaryArchive();
    ~BinaryArchive();

    int getVersion() const;
    std::string getHeader() const;

protected:
    int m_version;
    std::string m_header;

    void open(const char* filename,  ios_base::openmode mode);
    void open(const std::string &filename,  ios_base::openmode mode);

};

class IBinaryArchive : public BinaryArchive
{

public:
    IBinaryArchive(const std::string &filename);

    void open(const std::string &filename);

    std::istream & getStream();

    using BinaryArchive::operator>>;
    using BinaryArchive::read;

};

class OBinaryArchive : public BinaryArchive
{

public:
    OBinaryArchive(const std::string &filename, int version);
    OBinaryArchive(const std::string &filename, int version, const std::string &header);

    void open(const std::string &filename, int version);
    void open(const std::string &filename, int version, const std::string &header);

    std::ostream & getStream();

    using BinaryArchive::operator<<;
    using BinaryArchive::write;

};

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
    void read(std::istream &stream, const T &value);

    template<typename T>
    void read(std::istream &stream, const T &value, size_t size);

}

}

#endif

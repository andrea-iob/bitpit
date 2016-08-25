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

#include "binary_archive.hpp"

namespace bitpit {

/*!
    \ingroup Binary
    \class BinaryArchive
    \brief Base class for binary archives.

    The BinaryArchive is the blas class from which input and output binary
    archives are derived from.
*/

/*!
    Default constructor
*/
BinaryArchive::BinaryArchive()
    : m_version(VERSION_UNKNOW)
{
}

/*!
    Destructor
*/
BinaryArchive::~BinaryArchive()
{
    close();
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param filename is the filename to open
    \param mode specifies the opening mode
*/
void BinaryArchive::open(const char* filename,  ios_base::openmode mode)
{
    std::fstream::open(filename, std::ios::binary | mode);
    if (fail() && bad()) {
        return;
    }
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param filename is the filename to open
    \param mode specifies the opening mode
*/
void BinaryArchive::open(const std::string &filename,  ios_base::openmode mode)
{
    open(filename.c_str(), mode);
}

/*!
    Gets the version associated to the archive.

    \result The version associated to the archive.
*/
int BinaryArchive::getVersion() const
{
    return m_version;
}

/*!
    Gets the header associated to the archive.

    \result The header associated to the archive.
*/
std::string BinaryArchive::getHeader() const
{
    return m_header;
}

/*!
    \ingroup binary
    \class IBinaryArchive
    \brief Input binary archive.

    The IBinaryArchive class can read a binary archive, the binary archive
    has a fixed ASCII header.
*/

/*!
    Creates a new input archive.

    \param filename is the filename of the archive
*/
IBinaryArchive::IBinaryArchive(const std::string &filename)
{
    open(filename);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param filename is the filename of the archive
*/
void IBinaryArchive::open(const std::string &filename)
{
    // Reset the header
    m_header.clear();

    // Open the stream
    BinaryArchive::open(filename, std::ios::in | std::ios_base::binary);

    // Read the header
    char headerBuffer[HEADER_SIZE + 1];
    headerBuffer[HEADER_SIZE] = '\0';
    read(reinterpret_cast<char*>(&headerBuffer), HEADER_SIZE);
    m_header = std::string(headerBuffer);

    // Read the version
    read(reinterpret_cast<char*>(&m_version), sizeof(m_version));
}

/*!
    Get a reference to the input stream associated to the archive.

    \result A reference to the input stream associated to the archive.
*/
std::istream & IBinaryArchive::getStream()
{
    return *this;
}

/*!
    \ingroup binary
    \class OBinaryArchive
    \brief Output binary archive.

    The IBinaryArchive class can write a binary archive, the binary archive
    has a fixed ASCII header.
*/

/*!
    Creates a new output archive.

    \param filename is the filename of the archive
    \param version is the version of the archive
*/
OBinaryArchive::OBinaryArchive(const std::string &filename, int version)
{
    open(filename, version, "");
}

/*!
    Creates a new output archive.

    \param filename is the filename of the archive
    \param version is the version of the archive
    \param header is the header of the archive
*/
OBinaryArchive::OBinaryArchive(const std::string &filename, int version, const std::string &header)
{
    open(filename, version, header);
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param filename is the filename of the archive
    \param version is the version of the archive
*/
void OBinaryArchive::open(const std::string &filename, int version)
{
    open(filename, version, "");
}

/*!
    Opens the specified file, associating it with the stream object, so that
    input/output operations are performed on its content.

    \param filename is the filename of the archive
    \param version is the version of the archive
    \param header is the header of the archive, the length of the header is
    limited to HEADER_SIZE characters.
*/
void OBinaryArchive::open(const std::string &filename, int version, const std::string &header)
{
    // Open the stream
    BinaryArchive::open(filename, std::ios::out | std::ios_base::binary);

    // Write the header
    std::string archiveHeader(header);
    archiveHeader.resize(HEADER_SIZE, ' ');
    write(reinterpret_cast<const char*>(archiveHeader.data()), HEADER_SIZE);

    // Write the version
    write(reinterpret_cast<const char*>(&version), sizeof(version));
}

/*!
    Get a reference to the input stream associated to the archive.

    \result A reference to the input stream associated to the archive.
*/
std::ostream & OBinaryArchive::getStream()
{
    return *this;
}

}

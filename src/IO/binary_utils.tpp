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

#ifndef __BITPIT_BINARY_UTILS_TPP__
#define __BITPIT_BINARY_UTILS_TPP__

namespace bitpit {

/*!
    @ingroup Binary
    @{
*/

namespace binary {

/*!
    Write the given data to the specified stream in binary format.

    \param stream is the stream to write to
    \param value is the data to write
*/
template<typename T>
void write(std::ostream &stream, const T &value)
{
    write(stream, &value, sizeof(T));
}

/*!
    Write the given data to the specified stream in binary format.

    \param stream is the stream to write to
    \param value is the data to write
    \param size is the size, expressed in bytes, of the data to write
*/
template<typename T>
void write(std::ostream &stream, const T &value, size_t size)
{
    stream.write(&value, size);
}

/*!
    Write the given data to the specified stream in binary format.

    \param stream is the stream to write to
    \param value is the data to write
*/
template<typename T>
void write(std::ostream &stream, const T *value)
{
    write(stream, value, sizeof(T));
}

/*!
    Write the given data to the specified stream in binary format.

    \param stream is the stream to write to
    \param value is the data to write
    \param size is the size, expressed in bytes, of the data to write
*/
template<typename T>
void write(std::ostream &stream, const T *value, size_t size)
{
    stream.write(reinterpret_cast<const char*>(value), size);
}

/*!
    Read the given data to the specified stream in binary format.

    \param stream is the stream to read from
    \param[out] value on output it will contain the read value
*/
template<typename T>
void read(std::istream &stream, T &value)
{
    read(stream, &value, sizeof(T));
}

/*!
    Read the given data to the specified stream in binary format.

    \param stream is the stream to read from
    \param[out] value on output it will contain the read value
    \param size is the size, expressed in bytes, of the data to read
*/
template<typename T>
void read(std::istream &stream, T &value, size_t size)
{
    read(stream, &value, size);
}

/*!
    Read the given data to the specified stream in binary format.

    \param stream is the stream to read from
    \param[out] value on output it will contain the read value
*/
template<typename T>
void read(std::istream &stream, T *value)
{
    read(stream, value, sizeof(T));
}

/*!
    Read the given data to the specified stream in binary format.

    \param stream is the stream to read from
    \param[out] value on output it will contain the read value
    \param size is the size, expressed in bytes, of the data to read
*/
template<typename T>
void read(std::istream &stream, T *value, size_t size)
{
    stream.read(reinterpret_cast<char*>(value), size);
}

}

/*!
    @}
*/

}

#endif

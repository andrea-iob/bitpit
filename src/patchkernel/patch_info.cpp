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

#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif

#include "patch_info.hpp"
#include "patch_kernel.hpp"

namespace bitpit {

/*!
	\ingroup patchkernel
	\class PatchInfo

	\brief The PatchInfo class provides an interface for defining patch info.
*/

/*!
	Default constructor
*/
PatchInfo::PatchInfo()
{
}

/*!
	Destructor
*/
PatchInfo::~PatchInfo()
{
}

/*!
	Resets the information.
*/
void PatchInfo::reset()
{
	m_patch = nullptr;
	_reset();
}

/*!
	Extracts the information.

	\param patch is patch from which the informations will be extracted
*/
void PatchInfo::extract(PatchKernel const *patch)
{
	// Reset information
	reset();
	if (patch == nullptr) {
		return;
	}

	// Extract new information
	m_patch = patch;
	_extract(patch);
}

/*!
	Updates the information.

	\param patch is patch from which the informations will be extracted
*/
void PatchInfo::update()
{
	extract(m_patch);
}

/*!
	\ingroup patchkernel
	\class PatchTopologyInfo

	\brief Topology information about the patch.
*/

/*!
	Creates a new info.

	\param patch is patch from which the informations will be extracted
*/
PatchTopologyInfo::PatchTopologyInfo(PatchKernel const *patch)
{
	reset();

	m_global = false;

	extract(patch);
}

#if BITPIT_ENABLE_MPI==1
/*!
	Creates a new info.

	\param patch is patch from which the informations will be extracted
	\param global if set to true global information will be extracted
*/
PatchTopologyInfo::PatchTopologyInfo(PatchKernel const *patch, bool global)
{
	reset();

	m_global = global;

	extract(patch);
}
#endif

/*!
	Internal function to reset the information.
*/
void PatchTopologyInfo::_reset()
{
	m_nVertices       = 0;
	m_nOrphanVertices = 0;
	m_nFreeVertices   = 0;

	m_nFaces     = 0;
	m_nFreeFaces = 0;

	m_nCells       = 0;
	m_nOrphanCells = 0;
	m_nFreeCells   = 0;
}

/*!
	Internal function to extract the topology information from the patch.

	\param patch is patch from which the informations will be extracted
*/
void PatchTopologyInfo::_extract(PatchKernel const *patch)
{
	// Extract local information
	m_nVertices       = m_patch->getVertexCount();
	m_nOrphanVertices = m_patch->countOrphanVertices();
	m_nFreeVertices   = m_patch->countFreeVertices();

	m_nFaces     = m_patch->countFaces();
	m_nFreeFaces = m_patch->countFreeFaces();

	m_nCells       = m_patch->getCellCount();
	m_nOrphanCells = m_patch->countOrphanCells();
	m_nFreeCells   = m_patch->countFreeCells();

	if (!m_global || patch->getProcessorCount() <= 1) {
		return;
	}

	// Remove ghost data
	m_nCells -= m_patch->getGhostCount();
	for (auto cellItr = m_patch->ghostConstBegin(); cellItr != m_patch->ghostConstEnd(); ++cellItr) {
	}

	// Exchange data
	MPI_Comm communicator = m_patch->getCommunicator();

	long nLocalVertices = m_nVertices;
	MPI_Allreduce(&nLocalVertices, &m_nVertices, 1, MPI_LONG, MPI_SUM, communicator);

	long nOrphanVertices = m_nOrphanVertices;
	MPI_Allreduce(&nOrphanVertices, &m_nOrphanVertices, 1, MPI_LONG, MPI_SUM, communicator);

	long nFreeVertices = m_nFreeVertices;
	MPI_Allreduce(&nFreeVertices, &m_nFreeVertices, 1, MPI_LONG, MPI_SUM, communicator);

	long nLocalFaces = m_nFaces;
	MPI_Allreduce(&nLocalFaces, &m_nFaces, 1, MPI_LONG, MPI_SUM, communicator);

	long nLocalFreeFaces = m_nFreeFaces;
	MPI_Allreduce(&nLocalFreeFaces, &m_nFreeFaces, 1, MPI_LONG, MPI_SUM, communicator);

	long nLocalCells = m_nCells;
	MPI_Allreduce(&nLocalCells, &m_nCells, 1, MPI_LONG, MPI_SUM, communicator);

	long nLocalOrphanCells = m_nOrphanCells;
	MPI_Allreduce(&nLocalOrphanCells, &m_nOrphanCells, 1, MPI_LONG, MPI_SUM, communicator);

	long nLocalFreeCells = m_nFreeCells;
	MPI_Allreduce(&nLocalFreeCells, &m_nFreeCells, 1, MPI_LONG, MPI_SUM, communicator);
}

/*!
	Display information.

	\param[in,out] out output stream
	\param[in] padding (default = 0) number of leading spaces for
	formatted output
*/
void PatchTopologyInfo::display(std::ostream &out, unsigned int padding) const
{
	std::string indent = std::string(padding, ' ');

	// Vertex info
	out << indent<< "Vertices ..." << std::endl;
	out << indent<< "  - n. vertices .......... " << m_nVertices       << std::endl;
	out << indent<< "  - n. orphan vertices ... " << m_nOrphanVertices << std::endl;
	out << indent<< "  - n. free vertices ..... " << m_nFreeVertices   << std::endl;

	// Face info
	out << indent<< "Faces ..." << std::endl;
	out << indent<< "  - n. faces ............. " << m_nFaces          << std::endl;
	out << indent<< "  - n. free faces ........ " << m_nFreeCells      << std::endl;

	// Cell info
	out << indent<< "Cells ..." << std::endl;
	out << indent<< "  - n. cells ............. " << m_nCells          << std::endl;
	out << indent<< "  - n. orphan cells ...... " << m_nOrphanCells    << std::endl;
	out << indent<< "  - n. free cells ........ " << m_nFreeCells      << std::endl;
}

#if BITPIT_ENABLE_MPI==1
/*!
	\ingroup patchkernel
	\class PatchGlobalInfo

	\brief Global information about the patch.
*/

/*!
	Creates a new info.

	\param patch is patch from which the informations will be extracted
*/
PatchGlobalInfo::PatchGlobalInfo(PatchKernel const *patch)
{
	reset();

	extract(patch);
}

/*!
	Internal function to reset the information.
*/
void PatchGlobalInfo::_reset()
{
	m_cellLocalToGlobalMap.clear();
	m_nGlobalInternals.clear();
}

/*!
	Internal function to extract global information from the patch.

	\param patch is patch from which the informations will be extracted
*/
void PatchGlobalInfo::_extract(PatchKernel const *patch)
{
	long globalId;

	size_t exchangeDataSize = sizeof(globalId);
	std::unique_ptr<DataCommunicator> dataCommunicator;

	if (m_patch->getProcessorCount() > 1) {
		// Create the data communicator
		dataCommunicator = std::unique_ptr<DataCommunicator>(new DataCommunicator(m_patch->getCommunicator()));
		dataCommunicator->setTag(108);

		// Set and start the receives
		for (const auto entry : m_patch->getGhostExchangeTargets()) {
			const int rank = entry.first;
			const auto &list = entry.second;

			dataCommunicator->setRecv(rank, list.size() * exchangeDataSize);
			dataCommunicator->startRecv(rank);
		}
	}

	// Get the internal count of all the partitions
	m_nGlobalInternals.resize(m_patch->getProcessorCount());
	if (m_patch->getProcessorCount() > 1) {
		long nLocalInternals = m_patch->getInternalCount();
		MPI_Allgather(&nLocalInternals, 1, MPI_LONG, m_nGlobalInternals.data(), 1, MPI_LONG, m_patch->getCommunicator());
	} else {
		m_nGlobalInternals[0] = m_patch->getInternalCount();
	}

	// Evaluate the offset for the current partition
	long offset = 0;
	for (int i = 0; i < m_patch->getRank() - 1; ++i) {
		offset += m_nGlobalInternals[i];
	}

	// Evalaute the global id of the internal cells
	if (m_patch->getInternalCount() > 0) {
		auto cbeginInternals = m_patch->m_cells.cbegin();
		auto cendInternals   = ++m_patch->m_cells.getConstIterator(m_patch->m_lastInternalId);

		globalId = offset;
		for (auto cellItr = cbeginInternals; cellItr != cendInternals; ++cellItr) {
			m_cellLocalToGlobalMap.insert({cellItr.getId(), globalId++});
		}
	}

	// Communicate the global id of the ghost cells
	if (m_patch->getProcessorCount() > 1) {
		// Set and start the sends
		for (const auto entry : m_patch->getGhostExchangeSources()) {
			const int rank = entry.first;
			auto &list = entry.second;

			dataCommunicator->setSend(rank, list.size() * exchangeDataSize);
			SendBuffer &buffer = dataCommunicator->getSendBuffer(rank);
			for (long id : list) {
				buffer << m_cellLocalToGlobalMap.at(id);
			}
			dataCommunicator->startSend(rank);
		}

		// Receive the global ids of the ghosts
		int nCompletedRecvs = 0;
		while (nCompletedRecvs < dataCommunicator->getRecvCount()) {
			int rank = dataCommunicator->waitAnyRecv();
			const auto &list = m_patch->getGhostExchangeTargets(rank);

			RecvBuffer &buffer = dataCommunicator->getRecvBuffer(rank);
			for (long id : list) {
				buffer >> globalId;
				m_cellLocalToGlobalMap.insert({id, globalId});
			}

			++nCompletedRecvs;
		}

		// Wait for the sends to finish
		dataCommunicator->waitAllSends();
	}
}

/*!
	Gets the rank of the cell with the specified local id

	\param id is the local id of the cell
	\return The rank of the specified cell.
*/
int PatchGlobalInfo::getCellRankFromLocal(long id) const
{
	auto ghostOwnerItr = m_patch->m_ghostOwners.find(id);
	if (ghostOwnerItr != m_patch->m_ghostOwners.end()) {
		return ghostOwnerItr->second;
	} else {
		return m_patch->getRank();
	}
}

/*!
	Gets the rank of the cell with the specified global id

	\param id is the global id of the cell
	\return The rank of the specified cell.
*/
int PatchGlobalInfo::getCellRankFromGlobal(long id) const
{
	long offset = 0;
	for (int k = 0; k < m_patch->getProcessorCount(); ++k) {
		offset += m_nGlobalInternals[k];
		if (id < offset) {
			return k;
		}
	}

	return -1;
}

/*!
	Return the global id of the cell with the specified local id

	\param id is the local id of the cell
	\return The global id of the specified cell.
*/
long PatchGlobalInfo::getCellGlobalId(long id) const
{
	return m_cellLocalToGlobalMap.at(id);
}

/*!
	Return the map between local indexes and global indexes.

	\result The map between local indexes and global indexes.
*/
const std::unordered_map<long, long> & PatchGlobalInfo::getCellGlobalMap() const
{
	return m_cellLocalToGlobalMap;
}
#endif

}

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

#ifndef __BITPIT_PATCH_INFO_HPP__
#define __BITPIT_PATCH_INFO_HPP__

#include <unordered_map>
#include <vector>

namespace bitpit {

class PatchKernel;

class PatchInfo {

public:
	PatchInfo();

	virtual ~PatchInfo();

	void reset();
	void extract(PatchKernel const *patch);
	void update();

protected:
	PatchKernel const *m_patch;

	virtual void _reset() = 0;
	virtual void _extract(PatchKernel const *patch) = 0;

};

class PatchTopologyInfo : public PatchInfo {

public:
	PatchTopologyInfo(PatchKernel const *patch = nullptr);
#if BITPIT_ENABLE_MPI==1
	PatchTopologyInfo(PatchKernel const *patch, bool global);
#endif

	void display(std::ostream &out, unsigned int padding = 0) const;

protected:
	void _reset();
	void _extract(PatchKernel const *patch);

private:
	bool m_global;

	long m_nVertices;
	long m_nOrphanVertices;
	long m_nFreeVertices;

	long m_nFaces;
	long m_nFreeFaces;

	long m_nCells;
	long m_nOrphanCells;
	long m_nFreeCells;

};

#if BITPIT_ENABLE_MPI==1
class PatchGlobalInfo : public PatchInfo {

public:
	PatchGlobalInfo(PatchKernel const *patch = nullptr);

	int getCellRankFromLocal(long id) const;
	int getCellRankFromGlobal(long id) const;
	long getCellGlobalId(long id) const;
	const std::unordered_map<long, long> & getCellGlobalMap() const;

protected:
	void _reset();
	void _extract(PatchKernel const *patch);

private:
	std::unordered_map<long, long> m_cellLocalToGlobalMap;
	std::vector<long> m_nGlobalInternals;

};
#endif

}

#endif

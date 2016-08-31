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

#ifndef __BITPIT_PATCH_MANAGER_TPP__
#define __BITPIT_PATCH_MANAGER_TPP__

namespace bitpit {

/*!
	Create a new patch

	\tparam Args are the arguments to be passed to the patch constructor
	\param id is the type-id of the patch to create
	\param args are the arguments to be passed to the patch constructor
	\result A pointer to the newly created patch
*/
template<typename ...Args>
PatchKernel * PatchManager::create(int id, Args...args)
{
	return PatchFactory<Args...>::create(id, std::forward<Args>(args)...);
}

namespace patch {

	/*!
		Create a new patch

		\tparam Args are the arguments to be passed to the patch constructor
		\param id is the type-id of the patch to create
		\param args are the arguments to be passed to the patch constructor
		\result A pointer to the newly created patch
	*/
	template<typename ...Args>
	PatchKernel * create(int id, Args...args)
	{
		return manager().create(id, std::forward<Args>(args)...);
	}

}

}

#endif

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

#ifndef __BITPIT_PATCH_FACTORY_TPP__
#define __BITPIT_PATCH_FACTORY_TPP__

namespace bitpit {

template <typename ...Args>
BasePatchCreator<Args...>::~BasePatchCreator()
{
}

template <class Type, typename ...Args>
PatchCreator<Type, Args...>::PatchCreator(CreatorFunction creatorFunction)
	: m_creatorFunction(creatorFunction)
{
}

template <class Type, typename ...Args>
typename PatchCreator<Type, Args...>::PatchType* PatchCreator<Type, Args...>::create(Args...args) const
{
	if (m_creatorFunction) {
		return m_creatorFunction(std::forward<Args>(args)...);
	}

	return new Type(std::forward<Args>(args)...);
}


template<typename ...Args>
PatchFactory<Args...>::PatchFactory()
	: m_defaultCreator(0)
{
}

template<typename ...Args>
PatchFactory<Args...>::~PatchFactory()
{
	for (auto entry : m_creators) {
		delete entry.second;
	}
}

template<typename ...Args>
PatchFactory<Args...> & PatchFactory<Args...>::instance()
{
	static PatchFactory<Args...> factory;
	return factory;
}

template<typename ...Args>
PatchKernel * PatchFactory<Args...>::create(int id, Args... args)
{
	PatchFactory& factory = instance();
	std::cout << " create " << __PRETTY_FUNCTION__ << std::endl;
	std::cout << " create " << id << " :: " << factory.m_creators.size() << std::endl;

	if (!factory.isTypeRegistered(id)) {
		const BasePatchCreator<Args...> *defaultCreator = factory.getDefaultCreator();
		if (defaultCreator) {
			return defaultCreator->create(std::forward<Args>(args)...);
		} else {
			throw std::runtime_error ("Patch type has no creator");
		}
	}

	return factory.getCreator(id)->create(std::forward<Args>(args)...);
}

template<typename ...Args>
template<typename PatchType>
int PatchFactory<Args...>::registerDefaultType()
{
	return registerDefaultType(new PatchCreator<PatchType, Args...>());
}

template<typename ...Args>
int PatchFactory<Args...>::registerDefaultType(const BasePatchCreator<Args...> *creator)
{
	m_defaultCreator = creator;

	return 0;
}


template<typename ...Args>
template<typename PatchType>
int PatchFactory<Args...>::registerType(int id)
{
	return registerType(id, new PatchCreator<PatchType, Args...>());
}

template<typename ...Args>
int PatchFactory<Args...>::registerType(int id, const BasePatchCreator<Args...> *creator)
{
	std::cout << " register type " << __PRETTY_FUNCTION__ << std::endl;
	PatchFactory& factory = instance();

	std::cout << " register  Type" << id << std::endl;

	factory.unregisterType(id);
	m_creators.insert({id, creator});

	return 0;
}

template<typename ...Args>
int PatchFactory<Args...>::isTypeRegistered(int id)
{
	return (m_creators.count(id) != 0);
}

template<typename ...Args>
void PatchFactory<Args...>::unregisterType(int id)
{
	if (!isTypeRegistered(id)) {
		return;
	}

	delete m_creators.at(id);
	m_creators.erase(id);
}

template<typename ...Args>
const BasePatchCreator<Args...> * PatchFactory<Args...>::getDefaultCreator()
{
	return m_defaultCreator;
}

template<typename ...Args>
const BasePatchCreator<Args...> * PatchFactory<Args...>::getCreator(int id)
{
	return m_creators.at(id);
}

}

#endif

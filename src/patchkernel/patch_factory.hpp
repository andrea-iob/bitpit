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

#ifndef __BITPIT_PATCH_FACTORY_HPP__
#define __BITPIT_PATCH_FACTORY_HPP__

#include <unordered_map>

namespace bitpit {

class PatchKernel;

template <typename ...Args>
class BasePatchCreator {

public:
	virtual ~BasePatchCreator();

	virtual PatchKernel * create(Args...args) const = 0;

};

template <class Type, typename ...Args>
class PatchCreator : public BasePatchCreator<Args...> {

typedef Type* (*CreatorFunction) (Args...args);

public:
	typedef Type PatchType;

	PatchCreator(CreatorFunction creatorFunction = nullptr);

	virtual PatchType* create(Args...args) const;

private:
	CreatorFunction m_creatorFunction;

};

template <typename ...Args>
class PatchFactory {

private:
	PatchFactory();

	~PatchFactory();

	PatchFactory(const PatchFactory&);
	PatchFactory& operator=(const PatchFactory&);

public:
	static PatchFactory & instance();

	static PatchKernel * create(int id, Args... args);

	int registerDefaultType(const BasePatchCreator<Args...> *creator);
	int registerType(int id, const BasePatchCreator<Args...> *creator);
	void unregisterType(int id);
	int isTypeRegistered(int id);

	template<typename PatchType>
	int registerDefaultType();

	template<typename PatchType>
	int registerType(int id);

private:
	const BasePatchCreator<Args...> *m_defaultCreator;
	std::unordered_map<int, const BasePatchCreator<Args...> *> m_creators;

	const BasePatchCreator<Args...> * getDefaultCreator();
	const BasePatchCreator<Args...> * getCreator(int id);

};

}

#ifdef __BITPIT_ALLOW_PATCH_REGISTRATION__

#include "patch_factory.tpp"

#define REGISTER_PATCH_TYPE(id_type, id_args, PatchType, args...) \
int patch_factory_##PatchType_##id_args = PatchFactory<args>::instance().registerType(id_type, new PatchCreator<PatchType, ##args>()); \
template class PatchFactory<args>; \

#define UNREGISTER_PATCH_TYPE(id_type, args...) \
PatchFactory<args>::instance().unregisterType(id_type);

#else

#define REGISTER_PATCH_TYPE(id_type, id_args, PatchType, args...) \
extern template class PatchFactory<args>; \

#define UNREGISTER_PATCH_TYPE(id_type, args...) ;

#endif

#endif

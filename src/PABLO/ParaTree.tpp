/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
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

#ifndef __BITPIT_PABLO_PARA_TREE_TPP__
#define __BITPIT_PABLO_PARA_TREE_TPP__

namespace bitpit {

#if BITPIT_ENABLE_MPI==1
/** Communicate data provided by the user between the processes.
*/
template<class Impl>
void
ParaTree::communicate(DataCommInterface<Impl> & userData){
    //BUILD SEND BUFFERS
    DataCommunicator communicator(m_comm);
    size_t fixedDataSize = userData.fixedSize();
    std::map<int,u32vector >::iterator bitend = m_bordersPerProc.end();
    std::map<int,u32vector >::iterator bitbegin = m_bordersPerProc.begin();
    for(std::map<int,u32vector >::iterator bit = bitbegin; bit != bitend; ++bit){
        int  key = bit->first;
        const u32vector & pborders = bit->second;
        size_t buffSize = 0;
        size_t nofPbordersPerProc = pborders.size();
        if(fixedDataSize != 0){
            buffSize = fixedDataSize*nofPbordersPerProc;
        }
        else{
            for(size_t i = 0; i < nofPbordersPerProc; ++i){
                buffSize += userData.size(pborders[i]);
            }
        }
        //enlarge buffer to store number of pborders from this proc
        buffSize += sizeof(size_t);
        //build buffer for this proc
        communicator.setSend(key,buffSize);
        SendBuffer & sendBuffer = communicator.getSendBuffer(key);
        //store number of pborders from this proc at the begining
        sendBuffer << nofPbordersPerProc;

        //WRITE SEND BUFFERS
        for(size_t j = 0; j < nofPbordersPerProc; ++j){
            userData.gather(sendBuffer,pborders[j]);
        }
    }

    communicator.discoverRecvs();
    communicator.startAllRecvs();

    communicator.startAllSends();

    //READ RECEIVE BUFFERS
    int ghostOffset = 0;
    std::vector<int> recvRanks = communicator.getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());
    for(int rank : recvRanks){
        communicator.waitRecv(rank);
        RecvBuffer & recvBuffer = communicator.getRecvBuffer(rank);
        size_t nofGhostFromThisProc = 0;
        recvBuffer >> nofGhostFromThisProc;
        for(size_t k = 0; k < nofGhostFromThisProc; ++k){
            userData.scatter(recvBuffer, k+ghostOffset);
        }
        ghostOffset += nofGhostFromThisProc;
    }
    communicator.waitAllSends();
}

/** Distribute Load-Balancing the octants (with user defined weights) of the whole tree and data provided by the user
* over the processes of the job following the Morton order.
* Until loadBalance is not called for the first time the mesh is serial.
* Even distribute data provided by the user between the processes.
* \param[in] userData User interface to distribute the data during loadBalance.
* \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
*/
template<class Impl>
void
ParaTree::loadBalance(DataLBInterface<Impl> & userData, dvector* weight){
    loadBalance(&userData, weight);
}

/** Distribute Load-Balancing the octants (with user defined weights) of the whole tree and data provided by the user
* over the processes of the job following the Morton order.
* Until loadBalance is not called for the first time the mesh is serial.
* Even distribute data provided by the user between the processes.
* \param[in] userData User interface to distribute the data during loadBalance.
* \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
*/
template<class Impl>
void
ParaTree::loadBalance(DataLBInterface<Impl> * userData, dvector* weight){
    //Write info on log
    (*m_log) << "---------------------------------------------" << std::endl;
    (*m_log) << " LOAD BALANCE " << std::endl;

    m_lastOp = OP_LOADBALANCE;
    if (m_nproc>1){

        std::vector<uint32_t> partition(m_nproc);
        if (weight == NULL)
            computePartition(partition.data());
        else
            computePartition(partition.data(), weight);

        weight = NULL;

        privateLoadBalance(userData, partition.data());

        //Write info of final partition on log
        (*m_log) << " " << std::endl;
        (*m_log) << " Final Parallel partition : " << std::endl;
        (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
        for(int ii=1; ii<m_nproc; ii++){
            (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << std::endl;
        }
        (*m_log) << " " << std::endl;
        (*m_log) << "---------------------------------------------" << std::endl;

    }
    else{
        m_loadBalanceRanges.clear();

        (*m_log) << " " << std::endl;
        (*m_log) << " Serial partition : " << std::endl;
        (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
        (*m_log) << " " << std::endl;
        (*m_log) << "---------------------------------------------" << std::endl;
    }

}

/** Distribute Load-Balanced the octants (with user defined weights) of the whole tree and data provided by the user
* over the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
* The families of octants of a desired level are retained compact on the same process.
* Even distribute data provided by the user between the processes.
* \param[in] userData User interface to distribute the data during loadBalance.
* \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
* \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
*/
template<class Impl>
void
ParaTree::loadBalance(DataLBInterface<Impl> & userData, uint8_t & level, dvector* weight){
    loadBalance(&userData, level, weight);
}

/** Distribute Load-Balanced the octants (with user defined weights) of the whole tree and data provided by the user
* over the processes of the job. Until loadBalance is not called for the first time the mesh is serial.
* The families of octants of a desired level are retained compact on the same process.
* Even distribute data provided by the user between the processes.
* \param[in] userData User interface to distribute the data during loadBalance.
* \param[in] level Number of level over the max depth reached in the tree at which families of octants are fixed compact on the same process (level=0 is classic LoadBalance).
* \param[in] weight Pointer to a vector of weights of the local octants (weight=NULL is uniform distribution).
*/
template<class Impl>
void
ParaTree::loadBalance(DataLBInterface<Impl> * userData, uint8_t & level, dvector* weight){

    //Write info on log
    (*m_log) << "---------------------------------------------" << std::endl;
    (*m_log) << " LOAD BALANCE " << std::endl;

    m_lastOp = OP_LOADBALANCE;
    if (m_nproc>1){

        std::vector<uint32_t> partition(m_nproc);
        computePartition(partition.data(), level, weight);

        privateLoadBalance(userData, partition.data());

        //Write info of final partition on log
        (*m_log) << " " << std::endl;
        (*m_log) << " Final Parallel partition : " << std::endl;
        (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
        for(int ii=1; ii<m_nproc; ii++){
            (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << std::endl;
        }
        (*m_log) << " " << std::endl;
        (*m_log) << "---------------------------------------------" << std::endl;

    }
    else{
        m_loadBalanceRanges.clear();

        (*m_log) << " " << std::endl;
        (*m_log) << " Serial partition : " << std::endl;
        (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
        (*m_log) << " " << std::endl;
        (*m_log) << "---------------------------------------------" << std::endl;
    }

}

/**
* Distribute Load-Balancing octants and user data of the whole
* tree over the processes of the job following a given partition
* distribution. Until loadBalance is not called for the first time
* the mesh is serial.
* \param[in] userData User data that will be distributed among the
* processes.
* \param[in] partition Target distribution of octants over processes.
*/
template<class Impl>
void
ParaTree::privateLoadBalance(DataLBInterface<Impl> & userData, uint32_t* partition){
    privateLoadBalance(&userData, partition);
}

/**
* Distribute Load-Balancing octants and user data of the whole
* tree over the processes of the job following a given partition
* distribution. Until loadBalance is not called for the first time
* the mesh is serial.
* \param[in] userData User data that will be distributed among the
* processes.
* \param[in] partition Target distribution of octants over processes.
*/
template<class Impl>
void
ParaTree::privateLoadBalance(DataLBInterface<Impl> * userData, uint32_t* partition){

    // Update load balance ranges
    std::unordered_map<int, std::array<uint32_t, 2>> sendRanges = evalLoadBalanceSendRanges(partition);
    std::unordered_map<int, std::array<uint32_t, 2>> recvRanges = evalLoadBalanceRecvRanges(partition);

    m_loadBalanceRanges = LoadBalanceRanges(m_serial, sendRanges, recvRanges);

    // Load balance
    if(m_serial)
    {
        m_lastOp = OP_LOADBALANCE_FIRST;
        (*m_log) << " " << std::endl;
        (*m_log) << " Initial Serial distribution : " << std::endl;
        for(int ii=0; ii<m_nproc; ii++){
            (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]+1)) << std::endl;
        }

        uint32_t stride = 0;
        for(int i = 0; i < m_rank; ++i)
            stride += partition[i];
        LocalTree::octvector octantsCopy = m_octree.m_octants;
        LocalTree::octvector::const_iterator first = octantsCopy.begin() + stride;
        LocalTree::octvector::const_iterator last = first + partition[m_rank];

        m_octree.m_octants.assign(first, last);
        octvector(m_octree.m_octants).swap(m_octree.m_octants);
        m_octree.m_sizeOctants = m_octree.m_octants.size();

        first = octantsCopy.end();
        last = octantsCopy.end();

        if (userData) {
            userData->assign(stride,partition[m_rank]);
        }

        //Update and build ghosts here
        updateLoadBalance();
        computeGhostHalo();
    }
    else
    {
        (*m_log) << " " << std::endl;
        (*m_log) << " Initial Parallel partition : " << std::endl;
        (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(0))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[0]+1)) << std::endl;
        for(int ii=1; ii<m_nproc; ii++){
            (*m_log) << " Octants for proc	"+ std::to_string(static_cast<unsigned long long>(ii))+"	:	" + std::to_string(static_cast<unsigned long long>(m_partitionRangeGlobalIdx[ii]-m_partitionRangeGlobalIdx[ii-1])) << std::endl;
        }

        //empty ghosts
        m_octree.m_ghosts.clear();
        m_octree.m_sizeGhosts = 0;
        //compute new partition range globalidx
        assert(m_nproc > 0);
        std::vector<uint64_t> newPartitionRangeGlobalidx(m_nproc);
        for(int p = 0; p < m_nproc; ++p){
            newPartitionRangeGlobalidx[p] = 0;
            for(int pp = 0; pp <= p; ++pp)
                newPartitionRangeGlobalidx[p] += (uint64_t)partition[pp];
            --newPartitionRangeGlobalidx[p];
        }

        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        // Questa parte di codice potrebbe avere dei problemi di
        // overflow, a variabili int32_t vengono assegnati valori
        // uint32_t, inoltre si usa un int64_t per gli indice globali.
        // Il fatto di usare interi con segno è perchè si voule
        // gestire i casi in cui non c'è coda/testa (usando -1 e
        // nOctants + 1). Questo però dà problemi sia per i numeri
        // negativi sia per un possibile overflow nel calcolo di
        // (nOctants + 1).
        //
        // Io calcolerei solo headOffset e tailOffset (che sono sicuramente
        // uint32_t) e gestirei a parte i casi speciali headOffset == 0
        // e tailOffset == 0.
        //
        // Questo cidce commentato, assolutamente non testato, può essere
        // un punto di partenza.
//
//         // Compute head offset
//         uint32_t headOffset;
//         if (m_rank == 0) {
//             headOffset = 0;
//         } else if (newPartitionRangeGlobalidx[m_rank-1] == m_partitionRangeGlobalIdx[m_rank-1]) {
//             headOffset = 0;
//         } else {
//             headOffset = newPartitionRangeGlobalidx[m_rank-1] - m_partitionRangeGlobalIdx[m_rank-1];
//         }
//
//         // Compute tail ffset
//         uint32_t tailOffset;
//         if (m_rank == m_nproc - 1) {
//             tailOffset = 0;
//         } else if (m_rank == 0) {
//             if (getNumOctants() == newPartitionRangeGlobalidx[m_rank]) {
//                 tailOffset = 0;
//             } else {
//                 tailOffset = getNumOctants() - newPartitionRangeGlobalidx[m_rank] - 1;
//             }
//         else {
//             if (getNumOctants() > (newPartitionRangeGlobalidx[m_rank] - m_partitionRangeGlobalIdx[m_rank -1])) {
//                 tailOffset = getNumOctants() - newPartitionRangeGlobalidx[m_rank] - m_partitionRangeGlobalIdx[m_rank -1];
//             } else {
//                 tailOffset = 0;
//             }
//         }
//
//         //Initialize data communicator
//         DataCommunicator lbCommunicator(m_comm);
//
//         //Compute first predecessor and first successor to send buffers to
//         int firstPredecessor;
//         int firstSuccessor;
//         if (m_rank == 0) {
//             firstPredecessor = -1;
//             firstSuccessor   =  1;
//         } else {
//             uint64_t firstOctantGlobalIdx = m_partitionRangeGlobalIdx[m_rank-1] + 1;
//
//             firstPredecessor = -1;
//             if (headOffset > 0) {
//                 uint64_t globalLastHead = firstOctantGlobalIdx + headOffset - 1;
//                 for(int pre = m_rank - 1; pre >=0; --pre){
//                     if (globalLastHead <= newPartitionRangeGlobalidx[pre]) {
//                         firstPredecessor = pre;
//                     }
//                 }
//             }
//
//             firstSuccessor = m_nproc;
//             if (tailOffset > 0) {
//                 uint64_t globalFirstTail = firstOctantGlobalIdx + getNumOctants() - tailOffset;
//                 for(int post = m_rank + 1; post < m_nproc; ++post){
//                     if (globalFirstTail <= newPartitionRangeGlobalidx[post] && (uint64_t)globalFirstTail > newPartitionRangeGlobalidx[post-1]) {
//                         firstSuccessor = post;
//                     }
//                 }
//             }
//         }
//

        //find resident octants local offset lastHead(lh) and firstTail(ft)
        int32_t lh,ft;
        if(m_rank == 0)
            lh = -1;
        else{
            lh = (int32_t)(newPartitionRangeGlobalidx[m_rank-1] + 1 - m_partitionRangeGlobalIdx[m_rank-1] - 1 - 1);
        }
        if(lh < 0)
            lh = - 1;
        else if(lh > (int64_t) getNumOctants() - 1)
            lh = getNumOctants() - 1;

        if(m_rank == m_nproc - 1)
            ft = getNumOctants();
        else if(m_rank == 0)
            ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] + 1);
        else{
            ft = (int32_t)(newPartitionRangeGlobalidx[m_rank] - m_partitionRangeGlobalIdx[m_rank -1]);
        }
        if(ft > (int32_t)(getNumOctants() - 1))
            ft = getNumOctants();
        else if(ft < 0)
            ft = 0;

        //compute size Head and size Tail
        uint32_t headSize = (uint32_t)(lh + 1);
        uint32_t tailSize = (uint32_t)(getNumOctants() - ft);
        uint32_t headOffset = headSize;
        uint32_t tailOffset = tailSize;

        //Initialize data communicator
        DataCommunicator lbCommunicator(m_comm);

        //Compute first predecessor and first successor to send buffers to
        int64_t firstOctantGlobalIdx = 0;// offset to compute global index of each octant in every process
        int64_t globalLastHead = (int64_t) lh;
        int64_t globalFirstTail = (int64_t) ft; //lastHead and firstTail in global ordering
        int firstPredecessor = -1;
        int firstSuccessor = m_nproc;
        if(m_rank != 0){
            firstOctantGlobalIdx = (int64_t)(m_partitionRangeGlobalIdx[m_rank-1] + 1);
            globalLastHead = firstOctantGlobalIdx + (int64_t)lh;
            globalFirstTail = firstOctantGlobalIdx + (int64_t)ft;
            for(int pre = m_rank - 1; pre >=0; --pre){
                if((uint64_t)globalLastHead <= newPartitionRangeGlobalidx[pre])
                    firstPredecessor = pre;
            }
            for(int post = m_rank + 1; post < m_nproc; ++post){
                if((uint64_t)globalFirstTail <= newPartitionRangeGlobalidx[post] && (uint64_t)globalFirstTail > newPartitionRangeGlobalidx[post-1])
                    firstSuccessor = post;
            }
        }
        else if(m_rank == 0){
            firstSuccessor = 1;
        }

        int intBuffer = 0;
        int contatore = 0;
        //build send buffers from Head
        uint32_t nofElementsFromSuccessiveToPrevious = 0;
        if(headSize != 0){
            for(int p = firstPredecessor; p >= 0; --p){
                if(headSize < partition[p]){
                    intBuffer = (newPartitionRangeGlobalidx[p] - partition[p] );
                    intBuffer = abs(intBuffer);
                    nofElementsFromSuccessiveToPrevious = globalLastHead - intBuffer;
                    if(nofElementsFromSuccessiveToPrevious > headSize || contatore == 1)
                        nofElementsFromSuccessiveToPrevious  = headSize;

                    std::size_t buffSize = (std::size_t)nofElementsFromSuccessiveToPrevious * (std::size_t)Octant::getBinarySize();
                    //compute size of data in buffers
                    if (userData) {
                        if(userData->fixedSize()){
                            buffSize +=  userData->fixedSize() * nofElementsFromSuccessiveToPrevious;
                        }
                        else{
                            for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
                                buffSize += userData->size(i);
                            }
                        }
                    }
                    //add room for uint32_t, number of octants in this buffer
                    buffSize += sizeof(uint32_t);
                    lbCommunicator.setSend(p,buffSize);
                    SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                    //store the number of octants at the beginning of the buffer
                    sendBuffer << nofElementsFromSuccessiveToPrevious;

                    for(uint32_t i = (uint32_t)(lh - nofElementsFromSuccessiveToPrevious + 1); i <= (uint32_t)lh; ++i){
                        sendBuffer << m_octree.m_octants[i];
                        if (userData) {
                            userData->gather(sendBuffer,i);
                        }
                    }
                    if(nofElementsFromSuccessiveToPrevious == headSize)
                        break;

                    lh -= nofElementsFromSuccessiveToPrevious;
                    globalLastHead -= nofElementsFromSuccessiveToPrevious;
                    headSize = lh + 1;
                    ++contatore;
                }
                else{
                    nofElementsFromSuccessiveToPrevious = globalLastHead - (newPartitionRangeGlobalidx[p] - partition[p]);
                    std::size_t buffSize = (std::size_t)nofElementsFromSuccessiveToPrevious * (std::size_t)Octant::getBinarySize();
                    //compute size of data in buffers
                    if (userData) {
                        if(userData->fixedSize()){
                            buffSize +=  userData->fixedSize() * nofElementsFromSuccessiveToPrevious;
                        }
                        else{
                            for(int64_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
                                buffSize += userData->size(i);
                            }
                        }
                    }
                    //add room for uint32_t, number of octants in this buffer
                    buffSize += sizeof(uint32_t);
                    lbCommunicator.setSend(p,buffSize);
                    SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                    //store the number of octants at the beginning of the buffer
                    sendBuffer << nofElementsFromSuccessiveToPrevious;

                    for(int64_t i = lh - nofElementsFromSuccessiveToPrevious + 1; i <= lh; ++i){
                        sendBuffer << m_octree.m_octants[i];
                        if (userData) {
                            userData->gather(sendBuffer,i);
                        }
                    }
                    lh -= nofElementsFromSuccessiveToPrevious;
                    globalLastHead -= nofElementsFromSuccessiveToPrevious;
                    headSize = lh + 1;
                    if(headSize == 0)
                        break;
                }
            }

        }
        uint32_t nofElementsFromPreviousToSuccessive = 0;
        contatore = 0;
        //build send buffers from Tail
        if(tailSize != 0){
            for(int p = firstSuccessor; p < m_nproc; ++p){
                if(tailSize < partition[p]){
                    nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
                    if(nofElementsFromPreviousToSuccessive > tailSize || contatore == 1)
                        nofElementsFromPreviousToSuccessive = tailSize;

                    std::size_t buffSize = (std::size_t)nofElementsFromPreviousToSuccessive * (std::size_t)Octant::getBinarySize();
                    //compute size of data in buffers
                    if (userData) {
                        if(userData->fixedSize()){
                            buffSize +=  userData->fixedSize() * nofElementsFromPreviousToSuccessive;
                        }
                        else{
                            for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
                                buffSize += userData->size(i);
                            }
                        }
                    }
                    //add room for uint32_t, number of octants in this buffer
                    buffSize += sizeof(uint32_t);
                    lbCommunicator.setSend(p,buffSize);
                    SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                    //store the number of octants at the beginning of the buffer
                    sendBuffer << nofElementsFromPreviousToSuccessive;

                    for(uint32_t i = ft; i < ft + nofElementsFromPreviousToSuccessive; ++i){
                        sendBuffer << m_octree.m_octants[i];
                        if (userData) {
                            userData->gather(sendBuffer,i);
                        }
                    }
                    if(nofElementsFromPreviousToSuccessive == tailSize)
                        break;
                    ft += nofElementsFromPreviousToSuccessive;
                    globalFirstTail += nofElementsFromPreviousToSuccessive;
                    tailSize -= nofElementsFromPreviousToSuccessive;
                    ++contatore;
                }
                else{
                    nofElementsFromPreviousToSuccessive = newPartitionRangeGlobalidx[p] - globalFirstTail + 1;
                    uint32_t endOctants = ft + nofElementsFromPreviousToSuccessive - 1;
                    std::size_t buffSize = (std::size_t)nofElementsFromPreviousToSuccessive * (std::size_t)Octant::getBinarySize();
                    //compute size of data in buffers
                    if (userData) {
                        if(userData->fixedSize()){
                            buffSize +=  userData->fixedSize() * nofElementsFromPreviousToSuccessive;
                        }
                        else{
                            for(uint32_t i = ft; i <= endOctants; ++i){
                                buffSize += userData->size(i);
                            }
                        }
                    }
                    //add room for uint32_t, number of octants in this buffer
                    buffSize += sizeof(uint32_t);
                    lbCommunicator.setSend(p,buffSize);
                    SendBuffer &sendBuffer = lbCommunicator.getSendBuffer(p);
                    //store the number of octants at the beginning of the buffer
                    sendBuffer << nofElementsFromPreviousToSuccessive;

                    for(uint32_t i = ft; i <= endOctants; ++i ){
                        sendBuffer << m_octree.m_octants[i];
                        if (userData) {
                            userData->gather(sendBuffer,i);
                        }
                    }
                    ft += nofElementsFromPreviousToSuccessive;
                    globalFirstTail += nofElementsFromPreviousToSuccessive;
                    tailSize -= nofElementsFromPreviousToSuccessive;
                    if(tailSize == 0)
                        break;
                }
            }
        }

        lbCommunicator.discoverRecvs();
        lbCommunicator.startAllRecvs();
        lbCommunicator.startAllSends();

        uint32_t nofNewHead = 0;
        uint32_t nofNewTail = 0;

        //READ number of octants per sender
        std::vector<int> recvRanks = lbCommunicator.getRecvRanks();
        std::sort(recvRanks.begin(),recvRanks.end());
        std::vector<uint32_t> nofNewOverProcs(recvRanks.size());
        for(int rank : recvRanks){
            lbCommunicator.waitRecv(rank);
            RecvBuffer & recvBuffer = lbCommunicator.getRecvBuffer(rank);
            uint32_t nofNewPerProc;
            recvBuffer >> nofNewPerProc;
            nofNewOverProcs[rank] = nofNewPerProc;
            if(rank < m_rank)
                nofNewHead += nofNewPerProc;
            else if(rank > m_rank)
                nofNewTail += nofNewPerProc;
        }

        //MOVE RESIDENT TO BEGIN IN OCTANTS
        uint32_t resEnd = getNumOctants() - tailOffset;
        uint32_t nofResidents = resEnd - headOffset;
        uint32_t octCounter = 0;
        for(uint32_t i = headOffset; i < resEnd; ++i){
            m_octree.m_octants[octCounter] = m_octree.m_octants[i];
            if (userData) {
                userData->move(i,octCounter);
            }
            ++octCounter;
        }
        uint32_t newCounter = nofNewHead + nofNewTail + nofResidents;
        m_octree.m_octants.resize(newCounter, Octant(m_dim));
        m_octree.m_sizeOctants = m_octree.m_octants.size();
        if (userData) {
            userData->resize(newCounter);
        }
        //MOVE RESIDENTS IN RIGHT POSITION
        uint32_t resCounter = nofNewHead + nofResidents - 1;
        for(uint32_t k = 0; k < nofResidents ; ++k){
            m_octree.m_octants[resCounter - k] = m_octree.m_octants[nofResidents - k - 1];
            if (userData) {
                userData->move(nofResidents - k - 1,resCounter - k);
            }
        }

        //READ BUFFERS AND BUILD NEW OCTANTS
        newCounter = 0;
        bool jumpResident = false;
        for(int rank : recvRanks){
            RecvBuffer & recvBuffer = lbCommunicator.getRecvBuffer(rank);
            uint32_t nofNewPerProc = nofNewOverProcs[rank];
            if(rank > m_rank && !jumpResident){
                newCounter += nofResidents ;
                jumpResident = true;
            }
            for(int i = nofNewPerProc - 1; i >= 0; --i){
                recvBuffer >> m_octree.m_octants[newCounter];
                if (userData) {
                    userData->scatter(recvBuffer,newCounter);
                }
                ++newCounter;
            }
        }
        lbCommunicator.waitAllSends();
        octvector(m_octree.m_octants).swap(m_octree.m_octants);
        m_octree.m_sizeOctants = m_octree.m_octants.size();

        if (userData) {
            userData->shrink();
        }

        //Update and ghosts here
        updateLoadBalance();
        computeGhostHalo();
        uint32_t nofGhosts = getNumGhosts();
        if (userData) {
            userData->resizeGhost(nofGhosts);
        }

    }
}

#endif

}

#endif

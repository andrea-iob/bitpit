#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_PABLO.hpp"
#include <vector>

#include "DataCommInterface.hpp"

class UserDataLB : public bitpit::DataLBInterface<UserDataLB>{
public:
  static constexpr int N = 300;
  struct Data{
    int a;
    double b[N];
  };

  std::vector<Data>& data, ghostdata;

  size_t fixedSize() const
  {
    return sizeof(int) + N*sizeof(double);
  }
  //size_t size(const uint32_t e) const;
  void move(const uint32_t from, const uint32_t to)
  {
    data[to] = data[from];
  }

  template<class Buffer>
  void gather(Buffer & buff, const uint32_t e)
  {
    buff << data[e].a;
    for(int i=0; i<N; i++)
      buff << data[e].b[i];
  }

  template<class Buffer>
  void scatter(Buffer & buff, const uint32_t e)
  {
    buff >> data[e].a;
    for(int i=0; i<N; i++)
      buff >> data[e].b[i];
  }

  void assign(uint32_t stride, uint32_t length)
  {
    for(int i=0; i<length; i++)
      data[i] = data[i+stride];
  }
  void resize(uint32_t newSize)
  {
    data.resize(newSize);
  }
  void resizeGhost(uint32_t newSize)
  {
    ghostdata.resize(newSize);
  }
  void shrink() {}

  UserDataLB(std::vector<Data>& data_, std::vector<Data>& ghostdata_)
    : data(data_), ghostdata(ghostdata_)
  {}

  ~UserDataLB(){}
};

/**
 * Run the example.
 */
void run(int rank)
{
    /**<Instantation of a 3D pablo uniform object.*/
    bitpit::PabloUniform pablo11(3);
    /**<Set 2:1 balance for the octree.*/
    int idx = 0;
    pablo11.setBalance(idx,true);

    /** Set Periodic boundary conditions */
    pablo11.setPeriodic(0);
    pablo11.setPeriodic(2);
    pablo11.setPeriodic(4);

    /**<Compute the connectivity and write the octree.*/
    pablo11.computeConnectivity();

    /**<Refine globally one level and write the octree.*/
    pablo11.adaptGlobalRefine();
    pablo11.adaptGlobalRefine();
    pablo11.adaptGlobalRefine();

    pablo11.setMarker(10,1);
    pablo11.setMarker(15,-1);

    pablo11.adapt();

    assert( !pablo11.checkToAdapt() );

    uint8_t levels = 6;
    pablo11.loadBalance(levels);

    if( rank==0 )
    {
      pablo11.setMarker(10,1);
    }
    if( rank==2 )
    {
      pablo11.setMarker(15,-1);
    }
    pablo11.adapt();
    assert( !pablo11.checkToAdapt() );

    {
      std::vector<UserDataLB::Data> oct_data(pablo11.getNumOctants()), ghost_data(pablo11.getNumGhosts());
      UserDataLB data_lb(oct_data,ghost_data);
      pablo11.loadBalance(data_lb, levels);
    }
}

/*!
* Main program.
*/
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    int nProcs;
    int rank;
#if BITPIT_ENABLE_MPI==1
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    nProcs = 1;
    rank   = 0;
#endif

    // Initialize the logger
    bitpit::log::manager().initialize(bitpit::log::SEPARATE, false, nProcs, rank);
    bitpit::log::cout() << fileVerbosity(bitpit::log::NORMAL);
    bitpit::log::cout() << consoleVerbosity(bitpit::log::QUIET);

    // Run the example
    try {
        run(rank);
    } catch (const std::exception &exception) {
        bitpit::log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}

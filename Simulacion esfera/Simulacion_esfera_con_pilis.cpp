
#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <cmath>
#include <memory>
#include <stdlib.h>


//entonces toca revisar la conversion de unidades de longitud de los pili a las unidades de lattice, adicionalmente es necesario

using namespace plb;

typedef double T;

#define DESCRIPTOR descriptors::D3Q19Descriptor

std::string outDir("./tmp/");

struct SimulationParameters {
    
    
    std::vector<T> xDomain;                         // Extent in the x-direction of the physical simulation domain.
    std::vector<T> yDomain;                         // Extent in the y-direction of the physical simulation domain.
    std::vector<T> zDomain;                         // Extent in the z-direction of the physical simulation domain.
    
    std::vector<std::string> movingSurfaceFileNames;    // Files with the moving immersed surface geometries.
    plint numMovingSurfaces;                            // Number of moving immersed surfaces.
    plint numSurfaces;                                  // Number of all immersed surfaces.
    bool refineSurfaceMeshes;                           // Immersed surfaces are represented as triangle surface meshes.
    T targetMaxEdgeLength;                              // If surface mesh refinement has been selected, this parameter
    // defines the target value of the maximum edge length of the new
    // refined surface mesh.
    plint maxNumSurfaceRefinements;                     // This is the maximum number of surface refinements to be performed,
    // whether or not the target maximum triangle edge length has been
    // achieved.
    
    std::vector<ConnectedTriangleSet<T> > allSurfaces;
    
    
    T characteristicLength;                         // Length to define the resolution and the Reynolds number.
    plint resolution;                               // Total number of lattice nodes in the characteristic length.
    
    T dt;
    T masa;
    T momento_de_inercia;
    plint maxIter;                                  // Maximum number of iterations.
    plint statIter;                                 // Number of iterations for terminal output.
    plint outIter;                                  // Number of iterations for disk output.
    plint cpIter;                                   // Number of iterations for checkpointing.
    plint abIter;                                   // Number of iterations for checking for user-driven program abortion.
    plint startIter;                                // Number of initial iterations for smoothly increasing the inlet velocity.
    
    int ibIter;                                     // Iterations for the immersed boundary method.
    
    std::string abortFileName;                      // File for signaling program abortion.
    std::string xmlContinueFileName;                // XML file for restarting.
    std::string baseFileName;
    // Basename of the checkpoint files.
    
    T rho;                                          // Fluid density in physical units
    T nu;                                           // Fluid kinematic viscosity in physical units.
    T ambientPressure;                              // Absolute stagnation pressure in physical units.
    T cSmago;                                       // Smagorinsky parameter.
    Array<T,3> inletVelocity;                       // Inlet velocity vector in physical units.
    
    T radio_bacteria;
    T targetCSmago;
    bool useParallelIO;
    bool inicia_movimiento;
    
    
    
    Precision precision;                            // Precision for geometric operations.
    
    bool outputInDomain;                            // Save data on disk in a volume domain or not?
    Cuboid<T> outputCuboid;                         // Volume domain for disk output.

    
    T lx, ly, lz;                                   // Dimensions of the physical simulation domain.
    plint nx, ny, nz;                               // Grid dimensions of the simulation domain.
    T rho_LB;
    T ambientPressure_LB;
    
    Array<T,3> velocidad;
    Array<T,3> velocidad_angular;
    Array<T,3> eje_unitario;
    Array<T,3> vectores_unitarios;

    Array<T,3> inletVelocity_LB;
    std::vector<plint> startIds;
    std::vector<Array<T,3> > vertices;
    std::vector<T> areas;
    std::vector<int> flags;
    plint nextIter;
    T omega;                                        // Relaxation parameter for the fluid.
    T dx;                                           // Discrete spatial step.
    Array<T,3> physicalLocation;                    // Location of the physical domain.
    Array<T,3> centro_geometrico;
    // Location of the physical domain.

    Box3D lateral1, lateral2 ,lateral3, lateral4;
    Box3D interior;
    plint smallEnvelopeWidth;                       // Standard width.
    plint largeEnvelopeWidth;                       // For velocity because of immersed walls.
    bool saveDynamicContent;
    plint fileNamePadding;
    bool incompressibleModel;
    
    Box3D outputDomain;                             // Domains for disk output.
    
    //parametro para los pili
    
    T ka;
    T kAB;
    T kBA;
    T delta_X_AB;
    T delta_X_BA;
    T Ntot;
    T La0;
    T Lb0;
    T lp;
    T Kb_T;
    T L_inicial;
    
    
    int numero_de_pilis;
    int Na_inicial;
    int Nb_inicial;

};


T toLB(T physVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return (physVal - location[direction]) / dx;
}

Array<T,3> toLB(Array<T,3> const& physVal, T dx, Array<T,3> const& location)
{
    return (physVal - location) / dx;
}

void readUserDefinedSimulationParameters(std::string xmlInputFileName, SimulationParameters& param)
{
    XMLreader document(xmlInputFileName);
    
    document["geometry"]["simulationDomain"]["x"].read(param.xDomain);
    document["geometry"]["simulationDomain"]["y"].read(param.yDomain);
    document["geometry"]["simulationDomain"]["z"].read(param.zDomain);
    document["geometry"]["movingSurfaceFileNames"].read(param.movingSurfaceFileNames);
    param.numMovingSurfaces = param.movingSurfaceFileNames.size();
    param.numSurfaces = param.numMovingSurfaces;
    
    document["numerics"]["refineSurfaceMeshes"].read(param.refineSurfaceMeshes);
    if (param.refineSurfaceMeshes) {
        document["numerics"]["targetMaxEdgeLength"].read(param.targetMaxEdgeLength);
        document["numerics"]["maxNumSurfaceRefinements"].read(param.maxNumSurfaceRefinements);
        PLB_ASSERT(param.maxNumSurfaceRefinements > 0);
    }
    std::string precision;
   
        param.precision = DBL;
    
    
    document["numerics"]["masa"].read(param.masa);

    document["numerics"]["characteristicLength"].read(param.characteristicLength);
    document["numerics"]["resolution"].read(param.resolution);
    document["numerics"]["dt"].read(param.dt);

    document["numerics"]["maxIter"].read(param.maxIter);
    document["numerics"]["ibIter"].read(param.ibIter);
    
    document["numerics"]["startIter"].read(param.startIter);
    document["numerics"]["cSmago"].read(param.cSmago);
    document["numerics"]["inletVelocity"].read<T,3>(param.inletVelocity);
    
    
    document["numerics"]["abortFileName"].read(param.abortFileName);
    document["numerics"]["xmlContinueFileName"].read(param.xmlContinueFileName);
    document["numerics"]["baseFileName"].read(param.baseFileName);
    document["numerics"]["useParallelIO"].read(param.useParallelIO);
    
    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["nu"].read(param.nu);
    document["fluid"]["ambientPressure"].read(param.ambientPressure);
    
    document["output"]["statIter"].read(param.statIter);
    document["output"]["outIter"].read(param.outIter);
    document["output"]["cpIter"].read(param.cpIter);
    document["output"]["abIter"].read(param.abIter);
    
    
    
    
    
    document["output"]["outputInDomain"].read(param.outputInDomain);
    if (param.outputInDomain) {
        std::vector<T> x, y, z;
        document["output"]["outputDomain"]["x"].read(x);
        document["output"]["outputDomain"]["y"].read(y);
        document["output"]["outputDomain"]["z"].read(z);
        param.outputCuboid.lowerLeftCorner[0] = x[0];
        param.outputCuboid.lowerLeftCorner[1] = y[0];
        param.outputCuboid.lowerLeftCorner[2] = z[0];
        param.outputCuboid.upperRightCorner[0] = x[1];
        param.outputCuboid.upperRightCorner[1] = y[1];
        param.outputCuboid.upperRightCorner[2] = z[1];
    }
    
    
}

void computeOutputDomain(SimulationParameters& param)
{
    if (!param.outputInDomain) {
        return;
    }
    
    Array<T,3> llc = param.outputCuboid.lowerLeftCorner;
    Array<T,3> urc = param.outputCuboid.upperRightCorner;
    
    plint x0 = util::roundToInt(toLB(llc[0], 0, param.dx, param.physicalLocation));
    if (x0 < 0) {
        x0 = 0;
    } else if (x0 > param.nx-1) {
        x0 = param.nx-1;
    }
    plint y0 = util::roundToInt(toLB(llc[1], 1, param.dx, param.physicalLocation));
    if (y0 < 0) {
        y0 = 0;
    } else if (y0 > param.ny-1) {
        y0 = param.ny-1;
    }
    plint z0 = util::roundToInt(toLB(llc[2], 2, param.dx, param.physicalLocation));
    if (z0 < 0) {
        z0 = 0;
    } else if (z0 > param.nz-1) {
        z0 = param.nz-1;
    }
    
    plint x1 = util::roundToInt(toLB(urc[0], 0, param.dx, param.physicalLocation));
    if (x1 < 0) {
        x1 = 0;
    } else if (x1 > param.nx-1) {
        x1 = param.nx-1;
    }
    plint y1 = util::roundToInt(toLB(urc[1], 1, param.dx, param.physicalLocation));
    if (y1 < 0) {
        y1 = 0;
    } else if (y1 > param.ny-1) {
        y1 = param.ny-1;
    }
    plint z1 = util::roundToInt(toLB(urc[2], 2, param.dx, param.physicalLocation));
    if (z1 < 0) {
        z1 = 0;
    } else if (z1 > param.nz-1) {
        z1 = param.nz-1;
    }
    
    PLB_ASSERT(x1 >= x0 && y1 >= y0 && z1 >= z0);
    
    param.outputDomain = Box3D(x0, x1, y0, y1, z0, z1);
}


void defineOuterDomain(SimulationParameters& param)
{
    
    param.lateral1=  Box3D(0,      param.nx-1,       0     ,       0     ,         1,      param.nz-2);
    param.lateral2=  Box3D(0,      param.nx-1,   param.ny-1,   param.ny-1,         1,      param.nz-2);
    
    param.lateral3 = Box3D(0,      param.nx-1,       0     ,      param.ny-1,       0,             0    );
    param.lateral4 = Box3D(0,      param.nx-1,       0     ,      param.ny-1,   param.nz-1,   param.nz-1);
    
    param.interior= Box3D( 0,      param.nx-1,        1    ,      param.ny-2,      2  ,      param.nz-3);
    
    
}

void calculateDerivedSimulationParameters(SimulationParameters& param)
{

    param.Kb_T=4.10;
    param.Ntot=1000;
    param.Nb_inicial=7;
    param.Na_inicial=Ntot-param.Na_inicial;
    param.La0=0.75;
    param.Lb0=5.7;
    param.lp=2.7;
    param.ka=2.0;
    param.kAB=5.0e-2;
    param.kBA=700.0;
    param.delta_X_AB=0.5;
    param.delta_X_BA=-0.5;
    param.L_inicial=(param.Na_inicia)*param.La0+(param.Nb_inicial)*param.Lb0
    param.numero_de_pilis=10;
    
    param.inicia_movimiento=false;
    param.eje_unitario=Array<T,3>(0.,1.,0.);
    
    param.smallEnvelopeWidth = 1;
    param.largeEnvelopeWidth = 4;
    param.saveDynamicContent = true;
    param.fileNamePadding = 8;
    
    param.lx = param.xDomain[1] - param.xDomain[0];
    param.ly = param.yDomain[1] - param.yDomain[0];
    param.lz = param.zDomain[1] - param.zDomain[0];
    param.physicalLocation = Array<T,3>(param.xDomain[0], param.yDomain[0], param.zDomain[0]);
    
    param.dx = param.characteristicLength / (param.resolution - 1.0);
    
    param.nx = util::roundToInt(param.lx / param.dx) + 1;
    param.ny = util::roundToInt(param.ly / param.dx) + 1;
    param.nz = util::roundToInt(param.lz / param.dx) + 1;
    
    
    param.rho_LB = 1.0;
    param.ambientPressure_LB = (1.0/param.rho) * (param.dt*param.dt/(param.dx*param.dx)) * param.ambientPressure;
    param.inletVelocity_LB = param.inletVelocity * (param.dt / param.dx);
    
    
    param.velocidad=Array<T,3>(0.,0.,0.);
    param.velocidad_angular=Array<T,3>(0.0,0.0,0.0);
    param.vectores_unitarios=Array<T,3>(0.,0.,0.);


    
    if (param.refineSurfaceMeshes) {
        if (param.targetMaxEdgeLength < 0.0) {
            param.targetMaxEdgeLength = param.dx;
        }
    }
    

    T nu_LB = param.nu * param.dt / (param.dx * param.dx);
    //este parametro es el basicamente controla
    param.omega = 1.0 / (DESCRIPTOR<T>::invCs2 * nu_LB + 0.5);
    
    
    computeOutputDomain(param);
    defineOuterDomain(param);
    
    
    T NORMA_VELOCIDAD=sqrt(param.inletVelocity[0]*param.inletVelocity[0] +param.inletVelocity[1]*param.inletVelocity[1]+ param.inletVelocity[2]*param.inletVelocity[2]);
    
    
    pcout << "velocidad del sonido = " << nu_LB*(1.0/param.omega-0.5) << std::endl;
    pcout << "nu_LB = " << nu_LB<< std::endl;
    pcout << "numero de reynolds = " << NORMA_VELOCIDAD*param.characteristicLength/param.nu << std::endl;
    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "tiempo de relajacion = " << 1/param.omega<< std::endl;
    pcout << "velocidad de salida en unidades de lattice = " << param.inletVelocity_LB[0]<< std::endl;



}


void initializeImmersedSurfaceData(SimulationParameters& param)
{
    std::vector<std::string> allSurfaceFileNames;
    // Later we count on the fact that moving surfaces are included in the global vectors in the beginning.
    allSurfaceFileNames.insert(allSurfaceFileNames.end(), param.movingSurfaceFileNames.begin(), param.movingSurfaceFileNames.end());
    
    param.vertices.resize(0);
    param.areas.resize(0);
    param.flags.resize(0);
    param.startIds.resize(0);
    
    for (plint iSurface = 0; iSurface < param.numSurfaces; iSurface++) {
        TriangleSet<T> *surfaceTriangleSet = new TriangleSet<T>(allSurfaceFileNames[iSurface], param.precision);
        surfaceTriangleSet->writeBinarySTL(outDir + allSurfaceFileNames[iSurface]);
        
        
        if (param.refineSurfaceMeshes) {
            bool succeeded = surfaceTriangleSet->refineRecursively(param.targetMaxEdgeLength,
                                                                   param.maxNumSurfaceRefinements);
            if (!succeeded) {
                pcout << std::endl;
                pcout << "WARNING: The target maximum triangle edge length " << param.targetMaxEdgeLength
                << " for the immersed surface " << iSurface << std::endl
                << "         was not reached after " << param.maxNumSurfaceRefinements << " refinement iterations."
                << std::endl;
                pcout << std::endl;
            }
            FileName newSurfaceFileName(allSurfaceFileNames[iSurface]);
            newSurfaceFileName.setName(outDir + newSurfaceFileName.getName() + "_refined");
            surfaceTriangleSet->writeBinarySTL(newSurfaceFileName.get());
        }
        

        
        T maxEdgeLength = surfaceTriangleSet->getMaxEdgeLength();
        surfaceTriangleSet->scale(12.5/(100.0*param.dx));
        //380nm ->0.01
        Array <T,3> posicion_inicial=Array<T,3>(param.nx-1,param.ny-1,param.nz-1)*0.5;
        Cuboid<T> el_cubo = surfaceTriangleSet->getBoundingCuboid();
        Array<T,3> obstacleCenter = 0.5 * (el_cubo.lowerLeftCorner + el_cubo.upperRightCorner);
        param.centro_geometrico=posicion_inicial;
        surfaceTriangleSet->translate(posicion_inicial-obstacleCenter);
        
        
        
        param.radio_bacteria= 0.5 *std::abs(el_cubo.lowerLeftCorner[0] - el_cubo.upperRightCorner[0]);
        param.momento_de_inercia=param.masa*(2.0/5.0)*(param.radio_bacteria)*(param.radio_bacteria);
        

        pcout << "radio bacteria: " << param.radio_bacteria*param.dx<< std::endl;
        pcout << "momento de inercia: " << param.momento_de_inercia << std::endl;

        pcout << "locacion obstaculo  " << posicion_inicial[0]*param.dx  << "  "<< posicion_inicial[1]*param.dx  << " "<< posicion_inicial[2]*param.dx << std::endl;
        
        pcout << "centro geometrico:  " << obstacleCenter[0]*param.dx  << "  "<< obstacleCenter[1]*param.dx  << " "<< obstacleCenter[2]*param.dx << std::endl;


        
        ConnectedTriangleSet<T> connectedTriangleSet(*surfaceTriangleSet);
        delete surfaceTriangleSet;
        plint numVertices = connectedTriangleSet.getNumVertices();
        plint numTriangles = connectedTriangleSet.getNumTriangles();
        
        pcout << "The immersed surface " << iSurface <<" has " << numVertices
        << " vertices and " << numTriangles << " triangles." << std::endl;
        pcout << "The immersed surface " << iSurface <<" has a maximum triangle edge length of " << maxEdgeLength << std::endl;
        if (maxEdgeLength >= 4.0 * param.dx) {
            pcout << std::endl;
            pcout << "CAUTION: The maximum triangle edge length for the immersed surface " << iSurface << " is greater than "
            << " 4 times dx."
            << std::endl;
            pcout << "         The immersed boundary method will not work correctly. Surface refinement is necessary."
            << std::endl;
            pcout << std::endl;
            exit(1);
        } else if (maxEdgeLength > param.dx) {
            pcout << std::endl;
            pcout << "WARNING: The maximum triangle edge length for the immersed surface " << iSurface << " is greater than dx."
            << std::endl;
            pcout << "         The immersed boundary method might not work in an optimal way. Surface refinement is recommended."
            << std::endl;
            pcout << std::endl;
        }
        
        param.startIds.push_back(param.vertices.size());
        for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
            param.vertices.push_back(connectedTriangleSet.getVertex(iVertex));
            T area;
            Array<T,3> unitNormal;
            connectedTriangleSet.computeVertexAreaAndUnitNormal(iVertex, area, unitNormal);
            param.areas.push_back(area);
            param.flags.push_back(iSurface);
        }
        
        param.allSurfaces.push_back(connectedTriangleSet);  // Keep the initial surface topology.
    }
}


class VelocityFunction {
public:
    VelocityFunction(SimulationParameters& param_)
    : param(param_)
    { }
    
    Array<T,3> operator()(pluint id)
    {
        Array<T,3> const& position = param.vertices[id];
        Array<T,3> rotacion=getRotationalVelocity(position, param.velocidad_angular,param.vectores_unitarios, param.centro_geometrico);
        Array<T,3> translacion=param.velocidad;
        Array<T,3> velocidad_global=rotacion*param.radio_bacteria+translacion;
        return velocidad_global;

    }
    
private:
    SimulationParameters& param;
};

class pili
{
int Na, Nb;
Array<T,3> posicion_en_la_esfera;
Array<T,3> punto_de_anclaje;


    
public:
    void incializar(int,int);
    void se_anclo(Array<T,3>);
    void definir_posicion_sobre_la_esfera(Array<T,3>);
    Array<T,3> movimiento(Array<T,3>, Array<T,2>, SimulationParameters&);
    

    
};

void pili::inicializar(int a, int b)
{
    Na=a;
    Nb=b;
}

void pili::se_anclo(Array<T,3> d)
{
    punto_de_anclaje=d;

}


void pili::definir_posicion_sobre_la_esfera(Array<T,3> p)
{
    posicion_en_la_esfera=p;

}


Array<T,3> pili::movimiento(Array<T,3> desplazamiento, Array<T,2> p, SimulationParameters& param)
{   Array<T,3> posicion=punto_de_anclaje-posicion_en_la_esfera;
    posicion=posicion/sqrt(posicion[0]^2+posicion[1]^2+posicion[2]^2);
    L=desplazamiento[0]*posicion[0]+desplazamiento[0]*posicion[0]+desplazamiento[0]*posicion[0];

    double Lb=solucion_cubica(Na,Nb,L);
    double c=Lb0*Nb
    double fuerzas=fuerza(Lb,c)
    double La=nuevo_La(fuerzas,Na)

    pa=probabilidad_de_transicion(param.kAB,fuerzas,param.delta_X_AB)
    pb=probabilidad_de_transicion(param.kBA,fuerzas,param.delta_X_BA)
    
    
    if (pa>p[0] && (Na)>0):
        Na=Na-1;
        Nb=Nb+1;
        
    if (pb>p[1] && (Nb)>7) :
        Na=Na+1;
        Nb=Nb-1;
    
    return fuerzas*posicion;
    

}

class SurfaceNormalFunction {
public:
    SurfaceNormalFunction(SimulationParameters& param_)
    : param(param_)
    { }
    
    Array<T,3> operator()(pluint id)
    {
        plint iSurface = -1;
        for (plint i = 0; i < param.numSurfaces - 1; i++) {
            if ((plint) id >= param.startIds[i] && (plint) id < param.startIds[i+1]) {
                iSurface = i;
                break;
            }
        }
        if (iSurface == -1) {
            iSurface = param.numSurfaces - 1;
        }
        PLB_ASSERT(iSurface >= 0);
        
        plint localId = id - param.startIds[iSurface];
        T area;
        Array<T,3> unitNormal;
        param.allSurfaces[iSurface].computeVertexAreaAndUnitNormal(localId, area, unitNormal, &param.vertices,
                                                                   param.startIds[iSurface]);
        
        return unitNormal;
    }
    
private:
    SimulationParameters& param;
};


void saveMovingSurfaces(SimulationParameters& param, std::string baseName, plint iIter)
{
    plint numDigits = util::val2str(param.numMovingSurfaces).length();
    for (plint iMovingSurface = 0; iMovingSurface < param.numMovingSurfaces; iMovingSurface++) {
        TriangleSet<T>* triangleSet = param.allSurfaces[iMovingSurface].toTriangleSet(param.precision,
                                                                                      &param.vertices, param.startIds[iMovingSurface]);
        PLB_ASSERT(triangleSet != 0);
        triangleSet->scale(param.dx);
        triangleSet->translate(param.physicalLocation);
        
        std::string fname = createFileName(
                                           createFileName(baseName + "moving_surface_", iMovingSurface, numDigits+1)+"_", iIter, param.fileNamePadding)
        + ".stl";
        triangleSet->writeBinarySTL(fname);
        
        delete triangleSet;
    }
}

void readMovingSurfaces(SimulationParameters& param, std::string baseName, plint iIter)
{
    plint numDigits = util::val2str(param.numMovingSurfaces).length();
    for (plint iMovingSurface = 0; iMovingSurface < param.numMovingSurfaces; iMovingSurface++) {
        std::string fname = createFileName(
                                           createFileName(baseName + "moving_surface_", iMovingSurface, numDigits+1)+"_", iIter, param.fileNamePadding)
        + ".stl";
        
        TriangleSet<T> *surfaceTriangleSet = new TriangleSet<T>(fname, param.precision);
        surfaceTriangleSet->scale(12.5/(100.0*param.dx));
        surfaceTriangleSet->translate(-param.physicalLocation / param.dx);
        
        ConnectedTriangleSet<T> connectedTriangleSet(*surfaceTriangleSet);
        delete surfaceTriangleSet;
        plint numVertices = connectedTriangleSet.getNumVertices();
        
        plint startId = param.startIds[iMovingSurface];
        for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
            plint id = startId + iVertex;
            param.vertices[id] = connectedTriangleSet.getVertex(iVertex);
        }
    }
}

void createFluidBlocks(SimulationParameters& param, MultiBlockLattice3D<T,DESCRIPTOR>*& lattice,
                       MultiScalarField3D<T>*& rhoBar, MultiTensorField3D<T,3>*& j, MultiContainerBlock3D*& container,
                       std::vector<MultiBlock3D*>& lattice_rho_bar_j_arg)
{
    Dynamics<T,DESCRIPTOR> *dynamics = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);
    param.incompressibleModel = true;
    
    Box3D fullDomain(0, param.nx-1, 0, param.ny-1, 0, param.nz-1);
    lattice = generateMultiBlockLattice<T,DESCRIPTOR>(fullDomain, dynamics->clone(),
                                                      param.smallEnvelopeWidth).release();
    defineDynamics(*lattice, lattice->getBoundingBox(), dynamics->clone());
    delete dynamics;
    lattice->toggleInternalStatistics(false);
    
    rhoBar = generateMultiScalarField<T>((MultiBlock3D&) *lattice, param.largeEnvelopeWidth).release();
    rhoBar->toggleInternalStatistics(false);
    
    j = generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, param.largeEnvelopeWidth).release();
    j->toggleInternalStatistics(false);
    
    container = new MultiContainerBlock3D((MultiBlock3D&) *rhoBar);
    
    lattice_rho_bar_j_arg.clear();
    lattice_rho_bar_j_arg.push_back(lattice);
    lattice_rho_bar_j_arg.push_back(rhoBar);
    lattice_rho_bar_j_arg.push_back(j);
    integrateProcessingFunctional(
                                  new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(),
                                  lattice->getBoundingBox(), lattice_rho_bar_j_arg, 0);
    integrateProcessingFunctional(
                                  new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),lattice->getBoundingBox(), lattice_rho_bar_j_arg, 3); // Boundary conditions are executed at levels 1 and 2.
    
    
    
    
    
    // Integrate the immersed boundary processors in the lattice multi-block.
    
    std::vector<MultiBlock3D*> args;
    
    plint pl = 4;
    
    args.resize(0);
    args.push_back(container);
    integrateProcessingFunctional(
                                  new InstantiateImmersedWallDataWithIndexedTagging3D<T>(param.vertices, param.areas, param.flags),
                                  container->getBoundingBox(), *lattice, args, pl);
    pl++;
    
    for (plint i = 0; i < param.ibIter; i++) {
        args.resize(0);
        args.push_back(rhoBar);
        args.push_back(j);
        args.push_back(container);
        integrateProcessingFunctional(
                                      new IndexedInamuroIteration3D<T,VelocityFunction>(
                                                                                        VelocityFunction(param), 1.0 / param.omega, param.incompressibleModel),
                                      rhoBar->getBoundingBox(), *lattice, args, pl);
        pl++;
    }
    
}



void outerDomainBoundaryConditions(SimulationParameters const& param,
                                   MultiBlockLattice3D<T,DESCRIPTOR> *lattice,
                                   MultiScalarField3D<T> *rhoBar, MultiTensorField3D<T,3> *j,
                                   OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc)
{
    
    
    lattice->periodicity().toggleAll(true);
    rhoBar->periodicity().toggleAll(true);
    j->periodicity().toggleAll(true);
    
    
    bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral1, boundary::freeslip);
    bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral2, boundary::freeslip);
    bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral3);
    bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral4);
    
    
    
    
    
    
}

void initializeSimulation(SimulationParameters& param, bool continueSimulation, std::string xmlRestartFileName,
                          plint& iniIter, MultiBlockLattice3D<T,DESCRIPTOR> *lattice, std::vector<MultiBlock3D*>& lattice_rho_bar_j_arg,
                          std::vector<MultiBlock3D*>& checkpointBlocks)
{
    if (!continueSimulation) {
        applyProcessingFunctional(
                                  new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
                                  lattice->getBoundingBox(), lattice_rho_bar_j_arg);
    }
    pcout << std::endl;
    
    if (continueSimulation) {
        pcout << "Reading state of the simulation from file: " << xmlRestartFileName << std::endl;
        loadState(checkpointBlocks, iniIter, param.saveDynamicContent, xmlRestartFileName);
        lattice->getTimeCounter().resetTime(iniIter);
        readMovingSurfaces(param, param.baseFileName, iniIter);
        pcout << std::endl;
    }
    
    pcout << std::endl;
}

void writeVTK(SimulationParameters const& param, MultiBlockLattice3D<T,DESCRIPTOR> *lattice,
              plint iIter)
{
    T pressureScale = param.rho * (param.dx * param.dx) / (param.dt * param.dt) * DESCRIPTOR<T>::cs2;
    T pressureOffset = param.ambientPressure -
    param.rho_LB * param.rho * (param.dx * param.dx) / (param.dt * param.dt) * DESCRIPTOR<T>::cs2;
    
    if (param.outputInDomain) {
        std::string fname = createFileName(outDir + "domain_", iIter, param.fileNamePadding);
        VtkImageOutput3D<T> vtkOut(fname, param.dx, param.physicalLocation);
        
        std::auto_ptr<MultiScalarField3D<T> > p = computeDensity(*lattice, param.outputDomain);
        vtkOut.writeData<float>(*p, "pressure", pressureScale, pressureOffset);
        p.reset();
        
        std::auto_ptr<MultiTensorField3D<T,3> > v = computeVelocity(*lattice, param.outputDomain);
        vtkOut.writeData<float>(*computeNorm(*v), "velocityNorm", param.dx / param.dt);
        vtkOut.writeData<3,float>(*v, "velocity", param.dx / param.dt);
        
    }
    
    
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);

    
    
    // Command-line arguments
    
    if (argc != 2 && argc != 3) {
        pcout << "Usage: " << argv[0] << " xml-input-file-name [xml-continue-file-name]" << std::endl;
        exit(1);
    }
    
    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);
    
    std::string xmlRestartFileName;
    bool continueSimulation = false;
    if (argc == 3) {
        xmlRestartFileName = std::string(argv[2]);
        continueSimulation = true;
    }
    
    // Set the simulation parameters.
    
    SimulationParameters param;
    
    readUserDefinedSimulationParameters(xmlInputFileName, param);
    calculateDerivedSimulationParameters(param);
    
    global::IOpolicy().activateParallelIO(param.useParallelIO);


    
    std::cout.precision(10);
    std::scientific(std::cout);
    
    // Immersed surfaces.
    
    pcout << "Processing immersed surface geometries." << std::endl;
    
    initializeImmersedSurfaceData(param);
    
    // Fluid.
    
    pcout << "Generating fluid blocks." << std::endl;
    
    MultiBlockLattice3D<T,DESCRIPTOR> *lattice = 0;
    MultiScalarField3D<T> *rhoBar = 0;
    MultiTensorField3D<T,3> *j = 0;
    MultiContainerBlock3D *container = 0;
    std::vector<MultiBlock3D*> lattice_rho_bar_j_arg;
    
    createFluidBlocks(param, lattice, rhoBar, j, container, lattice_rho_bar_j_arg);
    
       // Boundary conditions.
    
    pcout << "Generating outer domain boundary conditions." << std::endl;
    
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    outerDomainBoundaryConditions(param, lattice, rhoBar, j, bc);
    delete bc;
    
    // Initialization.
    
    std::vector<MultiBlock3D*> checkpointBlocks;
    checkpointBlocks.push_back(lattice);
    checkpointBlocks.push_back(rhoBar);
    checkpointBlocks.push_back(j);
    
    plint iniIter = 0;
    
    initializeSimulation(param, continueSimulation, xmlRestartFileName, iniIter, lattice, lattice_rho_bar_j_arg,
                         checkpointBlocks);
    
    
    
    //generar pilis sobre la superficie de la esfera;
    
    pili Pilis[param.numero_de_pilis];
    double n;
    double vector;
    
    for (int j=0; j<param.numero_de_pilis; j=j+1)
    
    {   double x =drand48();
        double y =drand48();
        double z =drand48();
        n=x*x+y*y+z*z;
        vector=param.centro_geometrico-(Array<T,3>(x,y,z)*1/n);
        Pilis[j].definir_posicion_sobre_la_esfera(vector);
        Pilis[j].inicializar(param.Na_inicial,param.Nb_inicial);
        
    }
    
    
    
    
    // Starting iterations.
    pcout << "Starting simulation." << std::endl;
    bool stopExecution = false;
    
    
    initializeAtEquilibrium(*lattice, param.interior,param.rho_LB, Array<T,3>(0.,0.,0.));
    setBoundaryVelocity(*lattice, param.lateral4, Array<T,3>(0.,0.,0.));
    setBoundaryVelocity(*lattice, param.lateral3, param.inletVelocity_LB);
    
    
    param.nextIter = iniIter + 1;
    
    
    //la fuerza
    Array<T,3> force=Array<T,3>(0.,0.,0.);
    //el torque
    Array<T,3> torque=Array<T,3>(0.,0.,0.);
    
    
    
    
    Array<T,3> vectores_unitarios;
    for (plint iIter = iniIter; iIter < param.maxIter && !stopExecution; iIter++) {
        param.nextIter = iIter + 1;
        
        if (iIter % param.outIter == 0 || iIter == param.maxIter - 1) {
            writeVTK(param, lattice, iIter);
            saveMovingSurfaces(param, outDir, iIter);
            pcout << std::endl;}
        
        if ((param.cpIter > 0 && iIter % param.cpIter == 0 && iIter != iniIter) ||
            iIter == param.maxIter - 1) {
            pcout << "Saving the state of the simulation at iteration: " << iIter << std::endl;
            saveState(checkpointBlocks, iIter, param.saveDynamicContent, param.xmlContinueFileName,
                      param.baseFileName, param.fileNamePadding);
            saveMovingSurfaces(param, param.baseFileName, iIter);
            pcout << std::endl;
        }
        
        if (iIter % param.abIter == 0) {
            
            stopExecution = abortExecution(param.abortFileName, checkpointBlocks, iIter,
                                           param.saveDynamicContent, param.xmlContinueFileName,
                                           param.baseFileName, param.fileNamePadding);
            
            if (stopExecution) {
                saveMovingSurfaces(param, param.baseFileName, iIter);
                pcout << "Aborting execution at iteration: " << iIter << std::endl;
                pcout << std::endl;
            }
        }

    

        
        


        //aceleración en nuevo intervalo
        if (iIter>3000)
        {
            
            recomputeImmersedForce(SurfaceNormalFunction(param), param.omega, param.rho_LB, *lattice,*container, param.largeEnvelopeWidth, lattice->getBoundingBox(), param.incompressibleModel);
            force = -reduceImmersedForce<T>(*container, 0);
            torque = -reduceAxialTorqueImmersed<T>(*container,param.centro_geometrico, param.eje_unitario, 0);
            
        }
        
        
        
        //vector velocidad lineal después.
        param.velocidad=(force*1.0/param.masa)+param.velocidad;
        
        //vector velocidad angular después.
        param.velocidad_angular=(torque*1.0/param.momento_de_inercia)+param.velocidad_angular;
        
        
        T norma = sqrt(param.velocidad_angular[0]*param.velocidad_angular[0]+param.velocidad_angular[1]*param.velocidad_angular[1] +param.velocidad_angular[2]*param.velocidad_angular[2]);
        if (norma!=0.0)
        {
            param.vectores_unitarios=param.velocidad_angular/norma;
        }

        //mover el sistema en el instante
            
        param.centro_geometrico=param.centro_geometrico+param.velocidad;
        
        
        for (plint iMovingSurface = 0; iMovingSurface < param.numMovingSurfaces; iMovingSurface++)
        {
            
            plint numVertices = param.allSurfaces[iMovingSurface].getNumVertices();
            plint startId = param.startIds[iMovingSurface];
            
            for (plint iVertex = 0; iVertex < numVertices; iVertex++)
            {
                plint id = iVertex + startId;
                Array<T,3>& position = param.vertices[id];
                
                
                position = position+param.velocidad;
                if (norma!=0)
                {
                    
                    position = getRotatedPosition(position, param.velocidad_angular,param.vectores_unitarios, param.centro_geometrico);

                }
                

                
            }

            
        }
        
        
        pcout << iIter<<","<< torque[0] << "," << torque[1] << "," << torque[2]<<","<< force[0] << "," << force[1] << "," << force[2]<<","<< param.velocidad[0] << "," << param.velocidad[1] << "," << param.velocidad[2]<<","<<param.velocidad_angular[0] << "," <<  param.velocidad_angular[1] << "," << param.velocidad_angular[2]<<std::endl;
        


        

        

        global::timer("lb-iter").start();
        lattice->executeInternalProcessors(); // Execute all processors and communicate appropriately.
        global::timer("lb-iter").stop();
        lattice->incrementTime();
        
        
    }
    
    
    
    delete container;
    delete j;
    delete rhoBar;
    delete lattice;
    
    exit(0);
    
    
}

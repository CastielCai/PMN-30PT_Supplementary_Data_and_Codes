#include <iostream>
#include "mex.h"
//#include "omp.h"
#include <cmath>
#include "gpu/mxGPUArray.h"
#include <cuda_runtime_api.h>
#include <cufftXt.h>
#include <cuda_runtime.h>

#include "utility.cu"

#define M_PI 3.141592654f
#define eps 1.19209e-38f

// #include "helper_functions.h"
// #include "helper_cuda.h"
// #include "convolutionSeparable_common.h"


using namespace std;
#define NRANK 2

#if __CUDA_ARCH__ < 600
template <typename T>
__device__ double atomicAdd(T* address, T val)
{
    unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
        __double_as_longlong(val +
        __longlong_as_double(assumed)));

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);

    return __longlong_as_double(old);
}
#endif

__device__ static float atomicMax(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

__device__ static float atomicMin(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_i, assumed,
            __float_as_int(::fminf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

__device__ __forceinline__ float atomicMinFloat (float * addr, float value) {
        float old;
        old = (value >= 0) ? __int_as_float(atomicMin((int *)addr, __float_as_int(value))) :
             __uint_as_float(atomicMax((unsigned int *)addr, __float_as_uint(value)));

        return old;
}

__device__ __forceinline__ float atomicMaxFloat (float * addr, float value) {
    float old;
    old = (value >= 0) ? __int_as_float(atomicMax((int *)addr, __float_as_int(value))) :
         __uint_as_float(atomicMin((unsigned int *)addr, __float_as_uint(value)));

    return old;
}



template <typename T>
void __global__ update_positivity_l2(T*rec_new, int dimx, int dimy, T scale_factor){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    if (x<dimx ) {
        long long i = (long long)z*dimx*dimy + y*dimx +  x;
        rec_new[i] = max( 0.0, rec_new[i] ) * scale_factor;
    }
}

template <typename T>
void __global__ update_l2(T*rec_new, int dimx, int dimy, T scale_factor){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    if (x<dimx ) {
        long long i = (long long)z*dimx*dimy + y*dimx +  x;
        rec_new[i] = rec_new[i] * scale_factor;
    }
}

template <typename T>
void __global__ update_positivity(T*rec_new, int dimx, int dimy){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    if (x<dimx ) {
        long long i = (long long)z*dimx*dimy + y*dimx +  x;
        rec_new[i] = max( 0.0, rec_new[i] );
    }
}


template <typename T>
void __global__ compute_xyz_shift( const T*Matrix, const T* shift,  T*x_shift, T*y_shift, T*z_shift, int Num_pjs){
    int const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<Num_pjs ) {
        int index = 9*i;
        for (int j=0; j<4; j++){
            x_shift[4*i+j] = Matrix[index+0]*shift[2*j] + Matrix[index+3]*0.0 + Matrix[index+6]*shift[2*j+1] ;
            y_shift[4*i+j] = Matrix[index+1]*shift[2*j] + Matrix[index+4]*0.0 + Matrix[index+7]*shift[2*j+1] ;
            z_shift[4*i+j] = Matrix[index+2]*shift[2*j] + Matrix[index+5]*0.0 + Matrix[index+8]*shift[2*j+1] ;
        }
    }   
    
}

template <typename T>
void __global__ R1norm(const T *d_vec, double* R1, int N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) {        
        atomicAdd( R1 , (double)abs(d_vec[i]) );
    }
}

template <typename T>
void __global__ R1norm(const T *d_vec, T* R1, int N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) {        
        atomicAdd( R1 , abs(d_vec[i]) );
    }
}

template <typename T>
void __global__ computeCoords(T* coords, const int dimx, const int dimy, const int ncx, const int ncy, int ncz, long long N, long long starting_point){
    long long i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) { 
        i+=starting_point;
        coords[3*i]   =  int(           i%dimx ) - ncx  + 1 ;
        //coords[3*i+1] =  int(  ( i%(dimx*dimy) ) /dimx ) - ncy + 1;
        coords[3*i+1] =  ( int( i/dimx ) ) % dimy - ncy + 1;
        coords[3*i+2] =  int(    i/(dimx*dimy) ) - ncz + 1 ;
    }    
}
//xx_h = mod( ((0:692)'),7)-3;
//yy_h = floor( mod( (0:692)', 7*9) /7)-4; yy_h = mod( floor( (0:692)' /7), 9 )-4;
//zz_h = floor( (0:692)'/(7*9)) - 5;

template <typename T>
void __global__ setValue(T*residual, double val, int N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    //T o_ratio_inv = 1.0/o_ratio;
    if (i<N) {
        residual[i] = val;
    }
}

//static const int blockSize = 1024;
//static const int gridSize = 24; 



template <typename T>
void __global__ computeResidual(T*residual, const T* projections, const T scale, long long N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N) {
        residual[i] = residual[i]*scale - projections[i];
    }
}

template <typename T>
void __global__ compute_z_coord(T*z_coord,const T*Matrix, const int nrows, const int ncols, const int dimz,const T*nc, const int o_ratio,const T*x_s,const T*y_s,const T*z_s){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;    
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    
    //int origin_offset = 1;
    //long s ;  
    //#pragma omp parallel for default(shared) private(i,s) schedule(static)
    if (x<nrows){
        //const T coord_x = int(           i%nrows ) - nc[0]  + 1 ;
        //const T coord_y = ( int( i/nrows ) ) % ncols - nc[1] + 1;
        const T coord_x = x - nc[0] + 1;;
        const T coord_y = y - nc[1] + 1;
        const T coord_z = z - nc[2] + 1 ;

        //long index = i*3;
        //const T x_i = Matrix[0]*coord_x + Matrix[3]*coord_y + Matrix[6]*coord_z + nc[0];
        //const T y_i = Matrix[1]*coord_x + Matrix[4]*coord_y + Matrix[7]*coord_z + nc[1];
        const T z_c = Matrix[2]*coord_x + Matrix[5]*coord_y + Matrix[8]*coord_z;

        //T x_is = x_i - origin_offset;
        //T y_is = y_i - origin_offset;

        atomicMax( &z_coord[x + y*nrows]              , z_c );
        atomicMin( &z_coord[x + y*nrows + nrows*ncols], z_c );
        
    }
}

template <typename T>
void __global__ sum_slices( T* residual, const T* slices, int dimx, int dimy , int num_slices){
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;    
    if(x<dimx){
        double val=0;
        long long N = (long long int) dimx * dimy;
        for(int s=0; s<num_slices; s++){
            val += double( slices[x + y*dimx + s*N] );
        }
        residual[x + y*dimx] = val;
    }    
}

template <typename T>
void __global__ RadonTF(const T*data,const T*Matrix, const int dimx, const int dimy, const T*nc, const int dimx_ext, const int dimy_ext, const T* nc_ext, 
  const int o_ratio, const T*x_s,const T*y_s, T*slices, int k_min, int k_max, const T defocus_step){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;    
    int const x = blockIdx.x * blockDim.x + threadIdx.x;;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    
    long s ;  
    //#pragma omp parallel for default(shared) private(i,s) schedule(static)
    if (x<dimx_ext){
        const T & data_i = data[z*dimy_ext*dimx_ext + y*dimx_ext + x];
        //const T coord_x = int(           i%nrows ) - nc[0]  + 1 ;
        //const T coord_y = ( int( i/nrows ) ) % ncols - nc[1] + 1;
        const T coord_x = x - nc_ext[0] + 1;;
        const T coord_y = y - nc_ext[1] + 1;
        const T coord_z = z - nc_ext[2] + 1 ;

        //long index = i*3;
        const T x_i = Matrix[0]*coord_x + Matrix[3]*coord_y + Matrix[6]*coord_z + nc[0] - 1;
        const T y_i = Matrix[1]*coord_x + Matrix[4]*coord_y + Matrix[7]*coord_z + nc[1] - 1;
        const T z_c = Matrix[2]*coord_x + Matrix[5]*coord_y + Matrix[8]*coord_z;
        int k_slice = roundf( z_c/defocus_step );
        k_slice = max( min(k_slice, k_max), k_min);
        T * result = slices + (k_slice-k_min)* dimx*dimy ;

        for (s=0; s<o_ratio; s++){
            //for (i = 0; i < N; ++i) {
            T x_is = x_i + x_s[s];
            T y_is = y_i + y_s[s];

            // get coordinates of bounding grid locations
            long long x_1 = ( long long) floor(x_is) ;
            long long x_2 = x_1 + 1;
            long long y_1 = ( long long) floor(y_is) ;
            long long y_2 = y_1 + 1;
            
            if (x_1>=-1 && x_2<=dimx  &&  y_1>=-1 && y_2<=dimy ){ 
                T w_x1 = x_2 - x_is ;
                T w_x2 = 1   - w_x1;
                T w_y1 = y_2 - y_is ;
                T w_y2 = 1   - w_y1;            
                if (x_1==-1){
                    if (y_1==-1){
                        atomicAdd( &result[x_2 + y_2*dimx] , w_x2*w_y2 * data_i);
                    }
                    else if(y_2==dimy){
                        atomicAdd( &result[x_2 + y_1*dimx] , w_x2*w_y1 * data_i);
                    }
                    else{
                        atomicAdd( &result[x_2 + y_1*dimx] , w_x2*w_y1 * data_i);
                        atomicAdd( &result[x_2 + y_2*dimx] , w_x2*w_y2 * data_i);                    
                    }
                }
                else if (x_2==dimx){
                    if (y_1==-1){
                        atomicAdd( &result[x_1 + y_2*dimx] , w_x1*w_y2 * data_i);
                    }
                    else if(y_2==dimy){
                        atomicAdd( &result[x_1 + y_1*dimx] , w_x1*w_y1 * data_i);
                    }
                    else{
                        atomicAdd( &result[x_1 + y_1*dimx] , w_x1*w_y1 * data_i);
                        atomicAdd( &result[x_1 + y_2*dimx] , w_x1*w_y2 * data_i);                  
                    } 
                }
                else{
                    if (y_1==-1){
                        atomicAdd( &result[x_1 + y_2*dimx] , w_x1*w_y2 * data_i);
                        atomicAdd( &result[x_2 + y_2*dimx] , w_x2*w_y2 * data_i);
                    }
                    else if(y_2==dimy){
                        atomicAdd( &result[x_1 + y_1*dimx] , w_x1*w_y1 * data_i);
                        atomicAdd( &result[x_2 + y_1*dimx] , w_x2*w_y1 * data_i);
                    }
                    else{
                        atomicAdd( &result[x_1 + y_1*dimx] , w_x1*w_y1 * data_i);
                        atomicAdd( &result[x_1 + y_2*dimx] , w_x1*w_y2 * data_i);
                        atomicAdd( &result[x_2 + y_1*dimx] , w_x2*w_y1 * data_i);
                        atomicAdd( &result[x_2 + y_2*dimx] , w_x2*w_y2 * data_i);                  
                    }                               
                }
            }
        }
    }
}



template <typename T>
void __global__ RadonTpose_updateRec(const T* Matrix, const int dimx, const int dimy, const T* nc, const int dimx_ext, const int dimy_ext, const T* nc_ext, const T* slices, 
const int o_ratio, const T*x_s,const T*y_s, T* Rec, float dt, long long N, int k_min, int k_max, const T defocus_step){
    //long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    int const x = blockIdx.x * blockDim.x + threadIdx.x;
    int const y = blockIdx.y;
    int const z = blockIdx.z;
    int origin_offset = 1;
    long s;
    //#pragma omp parallel for default(shared) private(s) schedule(static)  
    if( x < dimx_ext ) {
        //const T coord_x = int(           i%nrows )   - nc[0] + 1 ;
        //const T coord_y = ( int( i/nrows ) ) % ncols - nc[1] + 1;
        //const T coord_z =  int(    i/(nrows*ncols) ) - nc[2] + 1 ;
        const T coord_x = x - nc_ext[0] + 1;;
        const T coord_y = y - nc_ext[1] + 1;
        const T coord_z = z - nc_ext[2] + 1 ;
        long long i = (long long)z*dimy_ext*dimx_ext + y*dimx_ext + x;

        //long index = i*3;
        const T x0  = Matrix[0]*coord_x + Matrix[3]*coord_y + Matrix[6]*coord_z + nc[0] - 1;
        const T y0  = Matrix[1]*coord_x + Matrix[4]*coord_y + Matrix[7]*coord_z + nc[1] - 1;
        const T z_c = Matrix[2]*coord_x + Matrix[5]*coord_y + Matrix[8]*coord_z;
        int k_slice = roundf( z_c/defocus_step );
        k_slice = max( min(k_slice, k_max), k_min);
        const T * data = slices + (k_slice-k_min)* dimx*dimy ;

        for (s=0; s<o_ratio; s++){
            T x_i = x0 + x_s[s];
            T y_i = y0 + y_s[s];
            // get coordinates of bounding grid locations
            long long x_1 = ( long long) floor(x_i) ;
            long long x_2 = x_1 + 1;
            long long y_1 = ( long long) floor(y_i) ;
            long long y_2 = y_1 + 1;
            
            // handle special case where x/y is the last element
            if ( x_i == dimx )   { x_2 -= 1; x_1 -= 1;}
            if ( y_i == dimx )   { y_2 -= 1; y_1 -= 1;}
            
            // return 0 for target values that are out of bounds
            /*
            if ( x_1>=0 && x_2<=(dimx-1)  && y_1>=0 && y_2<=(dimy - 1) ){
                // get the array values
                const T& f_11 = data[x_1 + y_1*dimx];
                const T& f_12 = data[x_1 + y_2*dimx];
                const T& f_21 = data[x_2 + y_1*dimx];
                const T& f_22 = data[x_2 + y_2*dimx];
                
                // compute weights
                T w_x1 = x_2 - x_i;
                T w_x2 = x_i - x_1;
                T w_y1 = y_2 - y_i;
                T w_y2 = y_i - y_1;
                
                T a,b;
                a = f_11 * w_x1 + f_21 * w_x2;
                b = f_12 * w_x1 + f_22 * w_x2;
                Rec[i] -= dt*(a * w_y1 + b * w_y2);
            }
            */

           if (x_1>=-1 && x_2<=dimx  &&  y_1>=-1 && y_2<=dimy ){             
                // compute weights
                T w_x1 = x_2 - x_i;
                T w_x2 = x_i - x_1;
                T w_y1 = y_2 - y_i;
                T w_y2 = y_i - y_1;

                // get the array values
                // const T f_11 = data[x_1 + y_1*dimx];
                // const T f_12 = data[x_1 + y_2*dimx];
                // const T f_21 = data[x_2 + y_1*dimx];
                // const T f_22 = data[x_2 + y_2*dimx];             
                if (x_1==-1){
                    if (y_1==-1){
                        Rec[i] -= dt * w_x2*w_y2* data[x_2 + y_2*dimx];
                    }
                    else if (y_2==dimy){
                        Rec[i] -= dt * w_x2*w_y1* data[x_2 + y_1*dimx];
                    }
                    else{
                        Rec[i] -= dt * ( w_x2*w_y2*data[x_2 + y_2*dimx] + w_x2*w_y1*data[x_2 + y_1*dimx] );
                    }
                }
                else if (x_2==dimx){
                    if (y_1==-1){
                        Rec[i] -= dt * w_x1*w_y2* data[x_1 + y_2*dimx];
                    }
                    else if (y_2==dimy){
                        Rec[i] -= dt * w_x1*w_y1* data[x_1 + y_1*dimx];
                    }
                    else{
                        Rec[i] -= dt * ( w_x1*w_y2*data[x_1 + y_2*dimx] + w_x1*w_y1*data[x_1 + y_1*dimx] );
                    }
                }
                else{
                    if (y_1==-1){
                        Rec[i] -= dt * ( w_x1*w_y2*data[x_1 + y_2*dimx] + w_x2*w_y2*data[x_2 + y_2*dimx] );
                    }
                    else if (y_2==dimy){
                        Rec[i] -= dt * ( w_x1*w_y1*data[x_1 + y_1*dimx] + w_x2*w_y1*data[x_2 + y_1*dimx] );
                    }
                    else{
                        Rec[i] -= dt * ( w_x1*w_y2*data[x_1 + y_2*dimx] + w_x2*w_y2*data[x_2 + y_2*dimx] +
                                         w_x1*w_y1*data[x_1 + y_1*dimx] + w_x2*w_y1*data[x_2 + y_1*dimx] );
                    }
                }
            }

        }
    }
}


template <typename T>
void __global__ avg_y(const T* d_slices, T*d_slices_tmp, const int dimx, const int dimy, const int num_slices){
    int i =  blockIdx.x * blockDim.x + threadIdx.x;
    int j =  blockIdx.y * blockDim.y + threadIdx.y;
    if(i<dimx && j<dimy){
        d_slices_tmp += (i+dimx*j);
        d_slices     += (i+dimx*j);
        if(j==0){
            for(int k=0; k<num_slices; k++){
                d_slices_tmp[0] = 0.8571429f* d_slices[0] + 0.1428571f* d_slices[dimx];
                d_slices_tmp += dimx*dimy;
                d_slices     += dimx*dimy;
            }      
        }
        else if(j==dimy-1){
            for(int k=0; k<num_slices; k++){
                d_slices_tmp[0] = 0.8571429f* d_slices[0] + 0.1428571f* d_slices[-dimx];
                d_slices_tmp += dimx*dimy;
                d_slices     += dimx*dimy;
            }        
        }
        else{
            for(int k=0; k<num_slices; k++){
                d_slices_tmp[0] =  0.75f* d_slices[0] + 0.125f* d_slices[dimx] + 0.125* d_slices[-dimx];
                d_slices_tmp += dimx*dimy;
                d_slices     += dimx*dimy;
            }        
        }
    }
}



template <typename T>
void __global__ compute_kR_rad2(T*kR_rad2, T dkx, const int nr, const int c_nr, const long long nr2){
    long long const ii = blockIdx.x * blockDim.x + threadIdx.x;
    if(ii<nr2){
        const int x = ii/nr - c_nr;
        const int y = ii%nr - c_nr;
        kR_rad2[ii] = (x*x + y*y)*dkx*dkx;
    }
}

template <typename T>
void __global__ compute_ape(bool*ape, const T*kR_rad2, T alpha2, const long long nr2){
    long long const ii = blockIdx.x * blockDim.x + threadIdx.x;
    if(ii<nr2){
        ape[ii] = kR_rad2[ii]<=alpha2;
    }
}

template <typename T>
void __global__ compute_probe_fun(float2* probe_fun, T* kR_rad2, bool* ape, T df_value, T lambda, long long N){
    long long ii =  blockIdx.x * blockDim.x + threadIdx.x;
    if(ii<N){
        float x=0,y=0;
        //float2 val = make_float2(1,0);
        if (ape[ii]) {
            float temp = kR_rad2[ii] * (df_value* 1e-6*M_PI/lambda);
            x = cosf(temp);
            y = sinf(temp);
        }
        probe_fun[ii].x = x;
        probe_fun[ii].y = y;
    }
}

template <typename T>
__global__ void cufftShift_2Da(T* data, int Nx, int Ny){
    // 2D Slice & 1D Line
    int Nxy = Nx * Ny;

    // Transformations Equations
    int sEq1 = Nx* ((Ny)/2) + (Nx)/2;//(Nxy + Nx) / 2;
    int sEq2 = sEq1 - Nx;

    // Thread Index (2D)
    int xIndex = blockIdx.x * blockDim.x + threadIdx.x;
    int yIndex = blockIdx.y * blockDim.y + threadIdx.y;

    // Thread Index Converted into 1D Index
    int index = yIndex*Nx + xIndex;

    T regTemp;

    if (xIndex < Nx/2) {
        if (yIndex < Ny / 2) {
            regTemp = data[index];
            // First Quad
            data[index] = data[index + sEq1];
            // Third Quad
            data[index + sEq1] = regTemp;
        }
    } 
    else {
        if (yIndex < Ny/2) {
            regTemp = data[index];
            // Second Quad
            data[index] = data[index + sEq2];
            // Fourth Quad
            data[index + sEq2] = regTemp;
        }
    }
}


template <typename T>
__global__ void fftshift_2D(const T* data_in, T*data_out, int Nx, int Ny, int num_slices){
    int ind_x = blockIdx.x * blockDim.x + threadIdx.x;
    int ind_y = blockIdx.y * blockDim.y + threadIdx.y;
    int ind_x_new = (ind_x + (Nx/2) ) % Nx;
    int ind_y_new = (ind_y + (Ny/2) ) % Ny;

    if (ind_x<Nx && ind_y<Ny) {
        for(int k=0; k<num_slices; k++){
            data_out[ ind_x_new + ind_y_new*Nx ] = data_in[ind_x + ind_y*Nx];
            data_in  += Nx*Ny;
            data_out += Nx*Ny;            
        }
    }
}


template <typename T>
__global__ void ifftshift_2D(const T* data_in, T*data_out, int Nx, int Ny, int num_slices){
    int ind_x = blockIdx.x * blockDim.x + threadIdx.x;
    int ind_y = blockIdx.y * blockDim.y + threadIdx.y;
    int ind_x_new = (ind_x + (Nx/2) ) % Nx;
    int ind_y_new = (ind_y + (Ny/2) ) % Ny;

    if (ind_x<Nx && ind_y<Ny) {
        for(int k=0; k<num_slices; k++){
            data_out[ind_x + ind_y*Nx] = data_in[ ind_x_new + ind_y_new*Nx ];            
            data_in  += Nx*Ny;
            data_out += Nx*Ny;
        }        
    }
}


void __global__ complexAbs2(float2 *d_in, float *d_out, long long N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) {        
        d_out[i] = d_in[i].x*d_in[i].x + d_in[i].y*d_in[i].y;
    }    
}

template <typename T>
void __global__ scalar_division(T *d_in, T scalar, long long N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) {        
        d_in[i] /= (scalar+eps);
    }
}

void __global__ Multiplication_half(const float2 *A, const float2 *B, float2*out, int dimx_half, int dimy){
    int ind_x = blockIdx.x * blockDim.x + threadIdx.x;
    int ind_y = blockIdx.y * blockDim.y + threadIdx.y;
    if (ind_x<dimx_half && ind_y<dimy ) {  
        long long i = dimx_half*ind_y + ind_x;      
        out[i].x = A[i].x*B[i].x - A[i].y*B[i].y;
        out[i].y = A[i].x*B[i].y + A[i].y*B[i].x;
        // if ( (ind_x+ind_y)%2){
        //     out[i].x = -out[i].x;
        //     out[i].y = -out[i].y;
        // }
    }    
}

void __global__ Multiplication(const float2 *A, const float2 *B, float2*out, long long N){
    long long const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<N ) {        
        out[i].x = A[i].x*B[i].x - A[i].y*B[i].y;
        out[i].y = A[i].x*B[i].y + A[i].y*B[i].x;
    }    
}


// (128,128) + (64,64) -> (192,192): 65+33=98, should be 65+32=97
// (128,128) + (65,65) -> (193,193): 65+33=98, should be 65+32=97
// (127,127) + (64,64) -> (191,191): 64+33=97, should be 64+32=96
// (127,127) + (65,65) -> (192,192): 64+33=97,

template <typename T>
void __host__ pad2D_zeros0(T*d_pjs, T*d_pjs_pad, int dimx, int dimy, int npad_x, int npad_y, int num_pjs ){
    int dimx_pad = dimx + npad_x;
    int dimy_pad = dimy + npad_y;
    int npad_x_half = npad_x/2;
    int npad_y_half = npad_y/2;
    if(npad_x%2 && dimx%2) npad_x_half++;
    if(npad_y%2 && dimy%2) npad_y_half++;
    
    for(int s=0; s<num_pjs; s++){
        for(int j=0; j<dimy; j++){
            cudaMemcpyAsync( d_pjs_pad+npad_x_half + (j+npad_y_half)*dimx_pad, d_pjs + j*dimx, 
                dimx*sizeof(T), cudaMemcpyDeviceToDevice);
        }
        d_pjs_pad  += dimx_pad*dimy_pad;
        d_pjs      += dimx*dimy;
    }
}


template <typename T>
void __host__ crop2D0(T*d_pjs, T*d_pjs_crop, int dimx_pad, int dimy_pad, int dimx, int dimy, int num_pjs ){
    int npad_x = dimx_pad - dimx;
    int npad_y = dimy_pad - dimy;
    int npad_x_half = npad_x/2;
    int npad_y_half = npad_y/2;
    if(npad_x%2 && dimx%2) npad_x_half++;
    if(npad_y%2 && dimy%2) npad_y_half++;

    for(int s=0; s<num_pjs; s++){
        for(int j=0; j<dimy; j++){
            cudaMemcpyAsync( d_pjs_crop + j*dimx, d_pjs + npad_x_half + (j+npad_y_half)*dimx_pad , 
                dimx*sizeof(T), cudaMemcpyDeviceToDevice);
        }
        d_pjs      += dimx_pad*dimy_pad;
        d_pjs_crop += dimx*dimy;
    }
}


int dimx;
int dimy;
int Num_pjs;

int dimx_pad;
int dimy_pad;
int dimx_pad_half;

int dimx_probe;
int dimx_probe_half;
int dimy_probe;


template <typename T>
void __global__ pad2D_zeros(T*d_pjs, T*d_pjs_pad, int dimx, int dimy, int dimx_pad, int dimy_pad, int npad_x_half, int npad_y_half, int num_pjs ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if(i<dimx && j<dimy){    
        d_pjs_pad += ( npad_x_half + npad_y_half*dimx_pad ); 
        for(int k=0; k<num_pjs; k++){
            d_pjs_pad[i + j*dimx_pad] = d_pjs[i + j*dimx];
            d_pjs_pad  += dimx_pad*dimy_pad;
            d_pjs      += dimx*dimy;
        }
    }
}

template <typename T>
void __global__ crop2D(T*d_pjs, T*d_pjs_crop, int dimx_pad, int dimy_pad, int dimx, int dimy, int npad_x_half, int npad_y_half, int num_pjs ){
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if(i<dimx && j<dimy){
        d_pjs += ( npad_x_half + npad_y_half*dimx_pad );  
        for(int k=0; k<num_pjs; k++){
            //d_pjs_crop[i + j*dimx] = d_pjs[npad_x_half + i + (j+npad_y_half)*dimx_pad];
            d_pjs_crop[i + j*dimx] = d_pjs[i + j*dimx_pad];
            d_pjs       += dimx_pad*dimy_pad;
            d_pjs_crop  += dimx*dimy;
        }
    }
}


template <typename T>
void __host__ FrequencyProbeFilter_forward(T* d_slices, T* kR_rad2, bool*ape, T* df_values, T lambda, T*hsize, int*dim_slices, int*dim_probe, cufftHandle& rfft_pad, cufftHandle& rifft_pad, cufftHandle& ifft_probe){    
    long long Nxy_probe     = dimx_probe*dimy_probe;
    long long Nxy_pad       = dimx_pad * dimy_pad;
    long long Nxy_pad_half  = dimx_pad_half*dimy_pad;

    // int& dimx         = dim_slices[0];
    // int& dimy         = dim_slices[1];
    int& num_slices   = dim_slices[2];
    //int dims[2]      = {dimy, dimx};
    //int dims_half[2] = {dimy, dimx/2+1};
    int dimx_half = dimx/2+1;
    long long Nxy      = ((long long int) dimx) * dimy;
    long long Nxy_half = (long long int) (dimx/2+1) * dimy;


    float2 *d_probe_1, *d_probe_2, *d_probe_3, *d_probe_4, *d_PSF, *d_FFT_A, *d_out;
    float *d_ker, *d_ker_pad, *d_slice_pad, *d_ker_pad_shift;
    cudaMalloc( &d_probe_1, Nxy_probe*sizeof(float2));
    cudaMalloc( &d_probe_2, Nxy_probe*sizeof(float2));
    cudaMalloc( &d_probe_3, Nxy_probe*sizeof(float2));
    cudaMalloc( &d_probe_4, Nxy_probe*sizeof(float2));
    cudaMalloc( &d_ker    , Nxy_probe*sizeof(float));

    cudaMalloc( &d_ker_pad        , Nxy_pad*sizeof(float));
    cudaMalloc( &d_ker_pad_shift  , Nxy_pad*sizeof(float));
    cudaMalloc( &d_slice_pad      , Nxy_pad*sizeof(float));
    cudaMemset( d_ker_pad   , 0.0f, Nxy_pad*sizeof(float));
    cudaMemset( d_slice_pad , 0.0f, Nxy_pad*sizeof(float));

    cudaMalloc( &d_PSF   , Nxy_pad_half*sizeof(float2));
    cudaMalloc( &d_FFT_A , Nxy_pad_half*sizeof(float2));
    cudaMalloc( &d_out   , Nxy_pad_half*sizeof(float2));

    // cufftHandle rfft_pad, ifft_probe, rifft_pad;
    // int dims_pad[2] = {dimy_pad, dimx_pad};
    // int dims_pad_half[2] = {dimy_pad, dimx_pad_half};
    // int dims_probe[2] = {dimy_probe, dimx_probe};
    // cufftPlanMany(&rfft_pad  , NRANK, dims_pad  , dims_pad     , 1, Nxy_pad     , dims_pad_half, 1, Nxy_pad_half, CUFFT_R2C, 1);
    // cufftPlanMany(&rifft_pad , NRANK, dims_pad  , dims_pad_half, 1, Nxy_pad_half, dims_pad     , 1, Nxy_pad     , CUFFT_C2R, 1);
    // cufftPlanMany(&ifft_probe, NRANK, dims_probe, dims_probe   , 1, Nxy_probe   , dims_probe   , 1, Nxy_probe   , CUFFT_C2C, 1);    
    
    for(int k=0; k<num_slices; k++){
        // abe_func     = kR_rad2* (-df_series(k) * 1e-6*pi/lambda);
        // probe_func   = exp(-1i*abe_func);
        // probe_func   = probe_func.*ape;
        compute_probe_fun<<< (Nxy_probe+255)/256, 256 >>>(d_probe_1, kR_rad2, ape, df_values[k], lambda, Nxy_probe);

        ifftshift_2D <<< dim3( (dimx_probe+15)/16, (dimy_probe+15)/16 ), dim3(16, 16) >>> (d_probe_1, d_probe_2, dimx_probe, dimy_probe, 1);
        
        if (cufftExecC2C(ifft_probe, d_probe_2, d_probe_3, CUFFT_INVERSE ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            return;	
        }

        fftshift_2D <<< dim3( (dimx_probe+15)/16, (dimy_probe+15)/16 ), dim3(16, 16) >>> (d_probe_3, d_probe_4, dimx_probe, dimy_probe, 1);  
        
        complexAbs2<<<(Nxy_probe+255)/256, 256>>>(d_probe_4, d_ker, Nxy_probe);        
        float norm_probe_value = sumGPUArray(d_ker, Nxy_probe);
        scalar_division <<< (Nxy_probe+255)/256, 256 >>> (d_ker, norm_probe_value, Nxy_probe);


        // PSF = fft2(probe): size change (nxn) -> (dimx_pad x dimy_pad)
        // d_ker will be destroyed
        pad2D_zeros0(d_ker, d_ker_pad, dimx_probe, dimy_probe, dimx-1, dimy-1, 1);

        ifftshift_2D <<< dim3( (dimx_pad+15)/16, (dimy_pad+15)/16 ), dim3(16, 16) >>> (d_ker_pad, d_ker_pad_shift, dimx_pad, dimy_pad, 1);
        if (cufftExecR2C( rfft_pad, d_ker_pad_shift, d_PSF ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            return;	
        }

        //float value1 = sumGPUArray(slices+k*Nxy, Nxy);

        // fft(probe): size (nxn) -> (dimx x dimy)
        // // B = fft2(A): slices will be destroyed
        pad2D_zeros0(d_slices+k*Nxy, d_slice_pad, dimx, dimy, dimx_probe-1, dimy_probe-1, 1);
        if (cufftExecR2C( rfft_pad, d_slice_pad, d_FFT_A ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            return;	
        }        

        // // F_A_blurry = fft2(A) .* fft2(probe_Kernel)
        //Multiplication_half <<< dim3( (dimx_pad_half+15)/16, (dimy_pad+15)/16 ) , dim3(16,16) >>> (d_FFT_A, d_PSF, d_out, dimx_pad_half, dimy_pad );
        Multiplication <<< (Nxy+255)/256, 256 >>> (d_FFT_A, d_PSF, d_out, Nxy_pad_half );

        // // A_blurry = ifft2( A .* PSF ): d_out will be destroyed        
        if (cufftExecC2R( rifft_pad, d_out, d_slice_pad ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            return;	
        }
        crop2D0(d_slice_pad, d_slices+k*Nxy, dimx_pad, dimy_pad, dimx, dimy, 1);
        
        float sc = (Nxy_pad+0.0f);
        scalar_division <<< (Nxy+255)/256, 256 >>> (d_slices+k*Nxy, sc, Nxy);

        //float value2 = sumGPUArray(slices+k*Nxy, Nxy);
        //cout << value1 << ", " << value2 << endl;
        cudaMemset( d_slice_pad , 0.0f, Nxy_pad*sizeof(float));
        cudaMemset( d_ker_pad   , 0.0f, Nxy_pad*sizeof(float));
    }
    
    // cufftDestroy(rfft_pad);
    // cufftDestroy(ifft_probe);
    // cufftDestroy(rifft_pad);
    
    cudaFree(d_probe_1);
    cudaFree(d_probe_2);
    cudaFree(d_probe_3);
    cudaFree(d_probe_4);
    cudaFree(d_ker);
    cudaFree(d_ker_pad);
    cudaFree(d_ker_pad_shift);
    cudaFree(d_slice_pad);
    cudaFree(d_out);
    cudaFree(d_PSF);
    cudaFree(d_FFT_A);
}



template <typename T>
void __host__ FrequencyProbeFilter(T* d_slices, T* kR_rad2, bool*ape, T* df_values, T lambda, T*hsize, int*dim_slices, int*dim_probe, cufftHandle& rfft_pad, cufftHandle& rifft_pad, cufftHandle& ifft_probe, float2* d_float2[7], float* d_float[4]){    
    long long Nxy_probe     = dimx_probe*dimy_probe;
    long long Nxy_pad       = dimx_pad * dimy_pad;
    long long Nxy_pad_half  = dimx_pad_half*dimy_pad;

    int& num_slices   = dim_slices[2];
    int dimx_half = dimx/2+1;
    long long Nxy      = ((long long int) dimx) * dimy;
    long long Nxy_half = (long long int) (dimx/2+1) * dimy;

    int npad_probe_x_half = (dimx-1)/2;
    int npad_probe_y_half = (dimy-1)/2;
    if( (dimx-1)%2 && dimx_probe%2) npad_probe_x_half++;
    if( (dimy-1)%2 && dimy_probe%2) npad_probe_y_half++;   

    int npad_x_half = (dimx_probe-1)/2;
    int npad_y_half = (dimy_probe-1)/2;
    if( (dimx_probe-1)%2 && dimx%2) npad_x_half++;
    if( (dimy_probe-1)%2 && dimy%2) npad_y_half++;      

    float2 *d_probe_1, *d_probe_2, *d_probe_3, *d_probe_4, *d_PSF, *d_FFT_A, *d_out;
    float *d_ker, *d_ker_pad, *d_slice_pad, *d_ker_pad_shift;
    d_probe_1 = d_float2[0];
    d_probe_2 = d_float2[1];
    d_probe_3 = d_float2[2];
    d_probe_4 = d_float2[3];

    d_PSF     = d_float2[4];
    d_FFT_A   = d_float2[5];
    d_out     = d_float2[6];
    
    d_ker           = d_float[0];
    d_ker_pad       = d_float[1];
    d_ker_pad_shift = d_float[2];
    d_slice_pad     = d_float[3];

    // cudaMalloc( &d_probe_1, Nxy_probe*sizeof(float2));
    // cudaMalloc( &d_probe_2, Nxy_probe*sizeof(float2));
    // cudaMalloc( &d_probe_3, Nxy_probe*sizeof(float2));
    // cudaMalloc( &d_probe_4, Nxy_probe*sizeof(float2));
    // cudaMalloc( &d_ker    , Nxy_probe*sizeof(float));
    // cudaMalloc( &d_ker_pad        , Nxy_pad*sizeof(float));
    // cudaMalloc( &d_ker_pad_shift  , Nxy_pad*sizeof(float));
    // cudaMalloc( &d_slice_pad      , Nxy_pad*sizeof(float));
    // cudaMalloc( &d_PSF   , Nxy_pad_half*sizeof(float2));
    // cudaMalloc( &d_FFT_A , Nxy_pad_half*sizeof(float2));
    // cudaMalloc( &d_out   , Nxy_pad_half*sizeof(float2));

    cudaMemset( d_ker_pad   , 0.0f, Nxy_pad*sizeof(float));
    cudaMemset( d_slice_pad , 0.0f, Nxy_pad*sizeof(float));

    for(int k=0; k<num_slices; k++){
        // abe_func     = kR_rad2* (-df_series(k) * 1e-6*pi/lambda);
        // probe_func   = exp(-1i*abe_func);
        // probe_func   = probe_func.*ape;
        compute_probe_fun<<< (Nxy_probe+255)/256, 256 >>>(d_probe_1, kR_rad2, ape, df_values[k], lambda, Nxy_probe);

        ifftshift_2D <<< dim3( (dimx_probe+15)/16, (dimy_probe+15)/16 ), dim3(16, 16) >>> (d_probe_1, d_probe_2, dimx_probe, dimy_probe, 1);
        
        if (cufftExecC2C(ifft_probe, d_probe_2, d_probe_3, CUFFT_INVERSE ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            cout <<"error: ifft_probe"<<endl;
            return;	
        }

        fftshift_2D <<< dim3( (dimx_probe+15)/16, (dimy_probe+15)/16 ), dim3(16, 16) >>> (d_probe_3, d_probe_4, dimx_probe, dimy_probe, 1);  
        
        complexAbs2<<<(Nxy_probe+255)/256, 256>>>(d_probe_4, d_ker, Nxy_probe);        
        float norm_probe_value = sumGPUArray(d_ker, Nxy_probe);
        scalar_division <<< (Nxy_probe+255)/256, 256 >>> (d_ker, norm_probe_value, Nxy_probe);


        // PSF = fft2(probe): size change (nxn) -> (dimx_pad x dimy_pad)
        // d_ker_pad_shift will be destroyed
        //pad2D_zeros(d_ker, d_ker_pad, dimx_probe, dimy_probe, dimx-1, dimy-1, 1);
        pad2D_zeros <<< dim3( (dimx_probe+15)/16, (dimy_probe+15)/16 ), dim3(16,16) >>> 
            (d_ker, d_ker_pad, dimx_probe, dimy_probe, dimx_pad, dimy_pad, npad_probe_x_half, npad_probe_y_half, 1);

        ifftshift_2D <<< dim3( (dimx_pad+15)/16, (dimy_pad+15)/16 ), dim3(16, 16) >>> (d_ker_pad, d_ker_pad_shift, dimx_pad, dimy_pad, 1);
        // error: need to fix
        if ( cufftExecR2C( rfft_pad, d_ker_pad_shift, d_PSF ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            cout <<"error: rfft_pad" <<endl;
            return;	
        }

        //float value1 = sumGPUArray(slices+k*Nxy, Nxy);

        // fft(probe): size (nxn) -> (dimx x dimy)
        // // B = fft2(A): d_slice_pad will be destroyed
        //pad2D_zeros(d_slices+k*Nxy, d_slice_pad, dimx, dimy, dimx_probe-1, dimy_probe-1, 1); 
        pad2D_zeros<<< dim3( (dimx+15)/16, (dimy+15)/16 ), dim3(16,16) >>>
            (d_slices+k*Nxy, d_slice_pad, dimx, dimy, dimx_pad, dimy_pad, npad_x_half, npad_y_half, 1);

        if (cufftExecR2C( rfft_pad, d_slice_pad, d_FFT_A ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            cout <<"error: rfft_pad"<<endl;
            return;	
        }        

        // // F_A_blurry = fft2(A) .* fft2(probe_Kernel)
        //Multiplication_half <<< dim3( (dimx_pad_half+15)/16, (dimy_pad+15)/16 ) , dim3(16,16) >>> (d_FFT_A, d_PSF, d_out, dimx_pad_half, dimy_pad );
        Multiplication <<< (Nxy_pad_half+255)/256, 256 >>> (d_FFT_A, d_PSF, d_out, Nxy_pad_half );

        // // A_blurry = ifft2( A .* PSF ): d_out will be destroyed        
        if (cufftExecC2R( rifft_pad, d_out, d_slice_pad ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
            fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
            cout <<"error: rifft_pad"<<endl;
            return;	
        }

        crop2D <<<dim3( (dimx+15)/16, (dimy+15)/16 ), dim3(16,16)>>>
            (d_slice_pad, d_slices+k*Nxy, dimx_pad, dimy_pad, dimx, dimy, npad_x_half, npad_y_half, 1);
        
        float sc = (Nxy_pad+0.0f);
        scalar_division <<< (Nxy+255)/256, 256 >>> (d_slices+k*Nxy, sc, Nxy);

        //float value2 = sumGPUArray(slices+k*Nxy, Nxy);
        //cout << value1 << ", " << value2 << endl;
        cudaMemset( d_slice_pad , 0.0f, Nxy_pad*sizeof(float));
        cudaMemset( d_ker_pad   , 0.0f, Nxy_pad*sizeof(float));
    }
    // cudaFree(d_probe_1a);
    // cudaFree(d_probe_2);
    // cudaFree(d_probe_3);
    // cudaFree(d_probe_4);
    // cudaFree(d_ker);
    // cudaFree(d_ker_pad);
    // cudaFree(d_ker_pad_shift);
    // cudaFree(d_slice_pad);
    // cudaFree(d_out);
    // cudaFree(d_PSF);
    // cudaFree(d_FFT_A);
}



void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) {
    int const threadsPerBlock = 256;
    int blocksPerGridPrj; // blocksPerGridRec;

    /* Initialize the MathWorks GPU API. */
    mxInitGPU();    
    cudaError_t err = cudaSuccess;
    err = cudaSetDevice(0);            // Set device 0 as current    
    if(err!=cudaSuccess){
        printf("cuda fail to set\n");
    }

    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    int device;
    for (device = 0; device < deviceCount; ++device) {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, device);
        printf("Device %d has compute capability %d.%d.\n",
        device, deviceProp.major, deviceProp.minor);
    }

    char const * const errId = "parallel:gpu:mexGPUExample:InvalidInput";
    char const * const errMsg = "Invalid input to MEX file.";    
    if ((nrhs!=7)&&(nrhs!=8) ){ //|| (nlhs!=1)  ) {
        cout << "error: number of parameter not matched"<<endl;
        mexErrMsgIdAndTxt(errId, errMsg); //!(mxIsGPUArray(prhs[0]))
    }
    /*
    * 0: projections
    * 1: rotmat
    * 2: dimz
    * 3: iterations
    * 4: step_size
    * 5: positivity?
    * 6: defocus_step
    * 7: Defocus_param
    * 8: semi_angle
    */

    // read inputs
    //1. projection
    const float  * projections  = mxGetSingles(prhs[0]);
    const size_t o_ratio=4;
    const mwSize* projsSize   = (mxGetDimensions(prhs[0]));
    dimx    = projsSize[0];
    dimy    = projsSize[1];
    Num_pjs = projsSize[2];

    //2. rotation matrix
    const float  * Matrix       = mxGetSingles(prhs[1]);
    
    //3. dimension of reconstruction
    const double * dim_extPtr = mxGetPr(prhs[2]);  
    const mwSize* dims_ext    = mxGetDimensions(prhs[2]);
    int num_dims_ext          = max(dims_ext[0], dims_ext[1]);
    if(num_dims_ext!=1 && num_dims_ext!=3){
        cout << "error: reconstruction dimension not matched"<<endl;
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    const int dimz            = num_dims_ext==1 ? dim_extPtr[0] :                dim_extPtr[2];
    const int dimx_ext        = num_dims_ext==1 ? dimx          : max(dimx, int( dim_extPtr[0] ) ); 
    const int dimy_ext        = num_dims_ext==1 ? dimy          : max(dimy, int( dim_extPtr[1] ) ); 
    const int dimz_ext        = dimz; 
    

    // 4. GD info
    const double * GD_info     = mxGetPr(prhs[3]);  
    const mwSize* dims_GD_info = mxGetDimensions(prhs[3]);
    int num_GD_info = max(dims_GD_info[0], dims_GD_info[1]);
    int iterations             = int   ( GD_info[0] ); 
    double step_size           = double( GD_info[1] ); 
    if(num_GD_info!=2 && num_GD_info!=Num_pjs+2){
        cout << "error: step size information not matched"<<endl;
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    const double dt =  (step_size/Num_pjs/dimz/o_ratio);
    float * theta_dt = new float[Num_pjs]{};
    if(num_GD_info==2){
        for(int i=0; i<Num_pjs; i++){
            theta_dt[i] = dt;
        }
    }
    else{
        const double* cos_theta    = GD_info + 2;
        for(int i=0; i<Num_pjs; i++){
            theta_dt[i] = dt*cos_theta[i];
        }
    }   

    // 5. constraints: positivity, is_y_avg, l2_reg
    const double * constraints     = mxGetPr(prhs[4]);  //mxGetLogicals(prhs[4]); 
    const mwSize* dims_constraints = mxGetDimensions(prhs[4]);
    int num_constraints   = max(dims_constraints[0], dims_constraints[1]);    
    if(num_constraints!=2 && num_constraints!=3){
        cout << "error: number of constraints not matched"<<endl;
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    const bool positivity = bool( constraints[0] ); 
    const bool is_y_avg   = constraints[1] != 0.0f;
    float l2_reg      = 0;
    if(num_constraints==3){
        l2_reg = constraints[2];
    }
    float scale_factor    = 1.0f / ( 1.0f + step_size*l2_reg );

    // const double * iterationPtr = mxGetPr(prhs[3]);  int iterations    = int( iterationPtr[0] ); 
    // const double * step_sizePtr = mxGetPr(prhs[4]);  double step_size  = double( step_sizePtr[0] ); 
    // const bool   * positivePtr  = mxGetLogicals(prhs[5]); bool positivity = bool( positivePtr[0] ); 
    
    // 6. defocus_info = [Voltage, pixelsize, nr, semi_angle, defocus_step];
    const double * defocus_info  = mxGetPr(prhs[5]);  
    double Voltage       = ( defocus_info[0] ); 
    double pixelsize     = ( defocus_info[1] ); 
    unsigned long nr     = (unsigned long) ( defocus_info[2] ); 
    double semi_angle    = ( defocus_info[3] ); 
    float defocus_step   = float( defocus_info[4] ); 
    float defocus_scale  = float( defocus_info[5] );

    // 7. defocus param
    const float  * defocus_param = mxGetSingles(prhs[6]);

    // print out params
    cout << "proj size  = (" << dimx << "," << dimy << "," << Num_pjs << ")" << endl;
    cout << "rec size   = (" << dimx_ext << "," << dimy_ext << "," << dimz_ext <<")" << endl;
    cout << "iterations = " << iterations << endl;
    cout << "step_size  = " << step_size << endl;
    cout << "using theta-step_size: " << (num_GD_info==Num_pjs+2) << std::endl;
    cout << "constraints:" <<endl;
    cout << "    positivity    : " << positivity << endl;
    cout << "    average on y  : " << is_y_avg << endl; 
    cout << "    l2 regularizer: " << l2_reg <<endl;


    const long long nrow_cols  = dimx*dimy;
    const long long nPjsPoints = dimx*dimy*Num_pjs;
    const long long recPoints = dimx_ext*dimy_ext*dimz_ext;
    const mwSize recSize[] = {dimx_ext,dimy_ext,dimz_ext};
    //cout << dimx <<", " << dimy <<", " << dimz <<", "  << endl;

    const mwSize* dims_Mat = (mxGetDimensions(prhs[1]));
    const mwSize Num_pjs2   = dims_Mat[2];
    const mwSize R1 = dims_Mat[0];
    const mwSize R2 = dims_Mat[1];

    const mwSize* dims_defocus_param = (mxGetDimensions(prhs[6]));
    const mwSize Num_pjs3   = dims_defocus_param[0];

    /*
    const mwSize* o_ratio_dim = mxGPUGetDimensions(xs2);
    const size_t o_ratio      = o_ratio_dim[0];
    const size_t Num_pjs2      = o_ratio_dim[1];                
    */
    if(Num_pjs2!=Num_pjs || Num_pjs3!=Num_pjs || R1!=3 || R2!=3 ){
        cout << "error: dimension not matched"<<endl;
        mexErrMsgIdAndTxt(errId, errMsg);
    }
    if(mxGetClassID(prhs[0]) == mxDOUBLE_CLASS || mxGetClassID(prhs[1]) == mxDOUBLE_CLASS ){     
        printf("can't work with double\n");
        return;
    }    

    // const mwSize ncx = int (floor(dimx/2.0)+1);
    // const mwSize ncy = int (floor(dimy/2.0)+1);
    // const mwSize ncz = int (floor(dimz/2.0)+1);

    
    // copy projections to GPU
    float * d_projections;
    cudaMalloc( &d_projections, nPjsPoints*sizeof(float) );
    cudaMemcpy(  d_projections, projections, nPjsPoints*sizeof(float), cudaMemcpyHostToDevice );

    // copy rotation matrix to GPU
    float * d_Matrix;
    cudaMalloc( &d_Matrix, 9*Num_pjs*sizeof(float) );
    cudaMemcpy( d_Matrix, Matrix, 9*Num_pjs*sizeof(float), cudaMemcpyHostToDevice );


    // plhs[0] is the returning reconstruction
    plhs[0] = mxCreateNumericArray(3, recSize, mxSINGLE_CLASS, mxREAL);
    float * rec   = (float*)mxGetSingles(plhs[0]);

    // plhs[1] is the returning calculated projection
    plhs[1] = mxCreateNumericArray(3, projsSize, mxSINGLE_CLASS, mxREAL);
    float* cal_proj   = (float*)mxGetSingles(plhs[1]);

    // create reconstruciton on GPU
    float * d_rec;
    cudaMalloc( &d_rec, recPoints*sizeof(float) );

    if(nrhs==8){
        const mwSize* recSize2   = mxGetDimensions(prhs[7]);
        const mwSize dimx2   = recSize2[0];
        const mwSize dimy2   = recSize2[1];
        const mwSize dimz2   = recSize2[2];
        if( dimx2!=dimx_ext || dimy2!=dimy_ext || dimz2!=dimz_ext ){
            cout << "error: reconstruction shape not matched"<<endl;
            mexErrMsgIdAndTxt(errId, errMsg);
        }
        float * rec_ori = (float*)mxGetSingles(prhs[7]);        
        //memcpy( h_rec, rec_ori, recPoints*sizeof(float) );
        cudaMemcpy( d_rec, rec_ori, recPoints*sizeof(float), cudaMemcpyHostToDevice);
    }     
    else{
        cudaMemset(d_rec, 0, recPoints*sizeof(float));
    }
    
    // create residual on GPU
    float * d_residual;
    cudaMalloc( &d_residual, nPjsPoints*sizeof(float) );
    /*
    mxGPUArray * residual = mxGPUCreateGPUArray(mwSize(3),
                        pjsSize, mxGPUGetClassID(projections),
                        mxREAL, MX_GPU_INITIALIZE_VALUES );
    */
    
    
    //blocksPerGrid    = (recPoints + threadsPerBlock - 1) / threadsPerBlock;
    blocksPerGridPrj = (nPjsPoints + threadsPerBlock - 1) / threadsPerBlock;
    //blocksPerGridRec = (recPoints + threadsPerBlock - 1) / threadsPerBlock;
    dim3 blocksPerGridRec2( (dimx_ext + threadsPerBlock - 1) / threadsPerBlock , dimy_ext, dimz_ext );

    //mxGPUArray const *support     = mxGPUCreateFromMxArray(prhs[7]);
    //bool const *d_support = (bool const *) (mxGPUGetDataReadOnly(support));            


    //float const *d_projections = (float const *) (mxGPUGetDataReadOnly(projections));
    //float const *d_Matrix      = (float const *) (mxGPUGetDataReadOnly(Matrix));
    //float const *d_Coord       = (float const *) (mxGPUGetDataReadOnly(Coord));
    //float const *d_xs         = (float const *) (mxGPUGetDataReadOnly(xs2));
    //float const *d_ys         = (float const *) (mxGPUGetDataReadOnly(ys2));
    //float *d_residual = (float *) (mxGPUGetData(residual));


    // compute rotated shift
    float shift[]  = {0.25,0.25, 0.25,-0.25, -0.25,0.25, -0.25,-0.25};
    float *shift_ptr;    
    cudaMalloc( (void**) &shift_ptr, 8*sizeof(float) );
    cudaMemcpy( shift_ptr, shift, 8*sizeof(float), cudaMemcpyHostToDevice);
    float * x_shift, *y_shift, * z_shift;
    cudaMalloc( (void**) &x_shift, 4*Num_pjs*sizeof(float) );
    cudaMalloc( (void**) &y_shift, 4*Num_pjs*sizeof(float) );
    cudaMalloc( (void**) &z_shift, 4*Num_pjs*sizeof(float) );
    compute_xyz_shift<<<2, threadsPerBlock>>>( d_Matrix, shift_ptr, x_shift, y_shift, z_shift, Num_pjs );
    float const *d_xs2         = (float  *) x_shift;
    float const *d_ys2         = (float  *) y_shift;
    float const *d_zs2         = (float  *) z_shift;

    // compute cartesian coordinates
    //computeCoords<<<blocksPerGridRec,threadsPerBlock>>>(d_Coord, dimx, dimy, ncx,ncy, ncz, recPoints, 0);

    // compute nc = [ncx,ncy,ncz]
    //const float nc_cpu[]  = { ncx,ncy,ncz}; 
    const float nc_ext_cpu[]  = { float(floor(dimx_ext/2.0)+1), float(floor(dimy_ext/2.0)+1), float(floor(dimz_ext/2.0)+1)}; 
    const float nc_cpu[]  = { float(floor(dimx/2.0)+1), float(floor(dimy/2.0)+1), float(floor(dimz/2.0)+1)};     
    float *d_nc, *d_nc_ext;
    cudaMalloc( &d_nc_ext, 3*sizeof(float) );
    cudaMalloc( &d_nc    , 3*sizeof(float) );
    cudaMemcpy( d_nc, nc_cpu        , 3*sizeof(float), cudaMemcpyHostToDevice ); 
    cudaMemcpy( d_nc_ext, nc_ext_cpu, 3*sizeof(float), cudaMemcpyHostToDevice ); 

    // compute norm of projection
    /*double pj_norm,  *pj_norm_pt;
    cudaMalloc( (void**)&pj_norm_pt, sizeof(double) );
    cudaMemset(pj_norm_pt, 0, 1*sizeof(double));   //setValue<<<1,1>>>(pj_pt, 0.0, 1);
    R1norm<<<blocksPerGridPrj, threadsPerBlock>>>(d_projections,  pj_norm_pt, nPjsPoints);
    cudaMemcpy( &pj_norm ,  pj_norm_pt, sizeof(double), cudaMemcpyDeviceToHost ) ;
    cudaFree( pj_norm_pt);
    */
    
    double pj_norm = sumArray(d_projections,nPjsPoints);

    // create temporary z_coord
    float * d_z_coord;
    cudaMalloc( &d_z_coord,       2*((long long) dimx)*dimy*sizeof(float) );
    cudaMemset(  d_z_coord, 0.0f, 2*((long long) dimx)*dimy*sizeof(float) );
    //cudaMemcpy(  d_projections, projections, nPjsPoints*sizeof(float), cudaMemcpyHostToDevice );
    int * pj_slices = new int[Num_pjs*2]{};
    int num_slices_max=0;
    for(int i=0; i<Num_pjs; i++){
        //cudaDeviceSynchronize();
        compute_z_coord<<<blocksPerGridRec2, threadsPerBlock>>>(d_z_coord, d_Matrix + i*9, dimx,dimy,dimz,d_nc,
            o_ratio, d_xs2+i*o_ratio,d_ys2+i*o_ratio,d_zs2+i*o_ratio);        
        double z_max = maxArray(d_z_coord , ( (long long) dimx )*dimy);
        double z_min = minArray(d_z_coord + dimx*dimy , ( (long long) dimx )*dimy);   
        int k_min = round(z_min/defocus_step);
        int k_max = round(z_max/defocus_step);
        num_slices_max = max(num_slices_max, k_max-k_min+1);
        pj_slices[2*i]   = k_min;
        pj_slices[2*i+1] = k_max;
        cudaMemset(d_z_coord, 0.0f, 2*dimx*dimy*sizeof(float) );
    }
    cudaFree(d_z_coord);    
    dim3 blocksPerGridSlice( (dimx + threadsPerBlock - 1) / threadsPerBlock , dimy );


    //defocus setting
    const double hh = 6.62e-34;
    const double ee = 1.6e-19;
    const double mm = 9.1e-31;
    const double cc = 3e8;
    const double lambda = 1e10*hh/sqrt(2*mm*ee*Voltage) / sqrt(1 + ee*Voltage/ (2.0*mm*cc*cc) );
    double kmax = 1.0e3/pixelsize*lambda/2;
    double dkx = 2*kmax/nr;

    // compute dimensions
    dimx_pad = dimx + nr;
    dimy_pad = dimy + nr;
    dimx_pad_half = dimx_pad/2+1;
    dimx_probe = nr+1;
    dimy_probe = nr+1;
    dimx_probe_half = dimx_probe/2+1;
    int dims[2]          = {dimy, dimx};
    int dims_pad[2]      = {dimy_pad, dimx_pad};
    int dims_probe[2]    = {dimy_probe, dimx_probe};
    int dims_half[2]     = {dimy, dimx/2+1};
    int dims_pad_half[2] = {dimy_pad, dimx_pad/2+1};
    long long Nxy = dimx*dimy;
    long long Nxy_pad = dimx_pad*dimy_pad;
    long long Nxy_probe = dimx_probe*dimy_probe;
    long long Nxy_half     = (long long) dimy*(dimx/2+1);
    long long Nxy_pad_half = (long long) dimy_pad*(dimx_pad/2+1);


    // compute kR_rad2
    float * d_kR_rad2;
    cudaMalloc( &d_kR_rad2,    Nxy_probe*sizeof(float) );
    cudaMemset(  d_kR_rad2, 0, Nxy_probe*sizeof(float) );
    compute_kR_rad2<<< (Nxy_probe+threadsPerBlock-1)/threadsPerBlock, threadsPerBlock >>> (d_kR_rad2, (float)dkx, dimx_probe, dimx_probe_half, Nxy_probe);


    // compute ape
    double alpha = semi_angle*1e3;
    float alpha2 = alpha*alpha;
    bool * d_ape;
    cudaMalloc( &d_ape,    Nxy_probe*sizeof(bool) );
    cudaMemset(  d_ape, 0, Nxy_probe*sizeof(bool) );
    compute_ape <<<(Nxy_probe+threadsPerBlock-1)/threadsPerBlock, threadsPerBlock >>> (d_ape, d_kR_rad2, alpha2, Nxy_probe);


    // test matlab probe function: plhs[2] = createProbeKernel
    mxArray *lhs[1], *rhs[5];
    double df_series = abs(2*defocus_step-0)*pixelsize;
    rhs[0] = mxCreateDoubleScalar(Voltage);
    rhs[1] = mxCreateDoubleScalar(pixelsize);
    rhs[2] = mxCreateDoubleScalar(nr);
    rhs[3] = mxCreateDoubleScalar(alpha);
    rhs[4] = mxCreateDoubleScalar(df_series);
    mexCallMATLAB(1, lhs, 5, rhs, "createProbeKernel");
    plhs[2] = lhs[0];


    //test fftshift: plhs[3] = fftshift_2D(projections)
    float *d_pj_in, *d_pj_out;
    cudaMalloc(&d_pj_in , dimx*dimy*sizeof(float));
    cudaMalloc(&d_pj_out, dimx*dimy*sizeof(float));
    cudaMemcpy(d_pj_in, d_projections, dimx*dimy*sizeof(float), cudaMemcpyHostToDevice);
    fftshift_2D<<< dim3( (dimx+15)/16, (dimy+15)/16 ), dim3(16, 16) >>>(d_pj_in, d_pj_out, dimx, dimy, 1);
    mwSize dim_1pj[] = {dimx, dimy};
    plhs[3] = mxCreateNumericArray(2, dim_1pj, mxSINGLE_CLASS, mxREAL);
    float *pj_out  = mxGetSingles(plhs[3]);
    cudaMemcpy(pj_out, d_pj_out, dimx*dimy*sizeof(float), cudaMemcpyDeviceToHost);    


    // call Matlab function myFrequencyProbeFilter
    // call FrequencyProbeFilter_forward(A, kR_rad2, ape, df_series, lambda , hsize)
    mxArray *LHS[1], *RHS[7];

    // copy kR_rad2
    const mwSize nrSize[] = {nr+1,nr+1};
    RHS[1] = mxCreateNumericArray(2, nrSize, mxSINGLE_CLASS, mxREAL);
    float * kR_rad2   = (float*)mxGetSingles(RHS[1]);
    cudaMemcpy( kR_rad2, d_kR_rad2, Nxy_probe*sizeof(float), cudaMemcpyDeviceToHost);

    // copy ape
    RHS[2] = mxCreateLogicalMatrix(nr+1, nr+1);
    bool * ape   = (bool*)mxGetLogicals(RHS[2]);
    cudaMemcpy( ape, d_ape, Nxy_probe*sizeof(bool), cudaMemcpyDeviceToHost);
    
    RHS[4] = mxCreateDoubleScalar(lambda);

    mwSize hsize_size[] = {1,2};
    RHS[5] = mxCreateNumericArray(2, hsize_size, mxSINGLE_CLASS, mxREAL);
    float * hsize   = (float*)mxGetSingles(RHS[5]);
    hsize[0] = nr;
    hsize[1] = nr;



    //test fft
    cufftHandle plan_fft2, plan_ifft2;
    if (cufftPlanMany(&plan_fft2, NRANK, dims, dims, 1, Nxy, dims_half, 1, Nxy_half, 
                    CUFFT_R2C, 1) != CUFFT_SUCCESS){
        fprintf(stderr, "CUFFT error: Plan creation failed");
        return;	
    }
    if (cufftPlanMany(&plan_ifft2, NRANK, dims, dims_half, 1, Nxy_half, dims, 1, Nxy, 
                    CUFFT_C2R, 1) != CUFFT_SUCCESS){
        fprintf(stderr, "CUFFT error: Plan creation failed");
        return;	
    }    
 

    float2 *d_pj_fft2;
    cudaMalloc(&d_pj_fft2 , (dimx/2+1)*dimy*sizeof(float2));
    if (cufftExecR2C( plan_fft2, d_pj_in, d_pj_fft2 ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
        fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
        return;	
    }
    mwSize dim_fft2[2] = {dimx/2+1, dimy};
    plhs[4] = mxCreateNumericArray(2, dim_fft2, mxSINGLE_CLASS, mxCOMPLEX);
    float2 *pj_fft  = (float2*)mxGetComplexSingles(plhs[4]);
    cudaMemcpy(pj_fft, d_pj_fft2, Nxy_half*sizeof(float2), cudaMemcpyDeviceToHost);  


    float *d_pj_ifft2;
    cudaMalloc(&d_pj_ifft2 , dimx*dimy*sizeof(float) );
    if (cufftExecC2R( plan_ifft2, d_pj_fft2, d_pj_ifft2 ) != CUFFT_SUCCESS){   //CUFFT_FORWARD
        fprintf(stderr, "CUFFT error: ExecC2C Inverse failed");
        return;	
    }
    mwSize dim_ifft2[2] = {dimx, dimy};
    plhs[5] = mxCreateNumericArray(2, dim_ifft2, mxSINGLE_CLASS, mxREAL);
    float *pj_ifft  = (float*)mxGetSingles(plhs[5]);
    cudaMemcpy(pj_ifft, d_pj_ifft2, Nxy*sizeof(float), cudaMemcpyDeviceToHost);      

    cudaFree(d_pj_in);
    cudaFree(d_pj_out);
    cudaFree(d_pj_fft2);
    cufftDestroy(plan_fft2);
    cufftDestroy(plan_ifft2);




    //fft plans
    cufftHandle rfft_pad, ifft_probe, rifft_pad;
    cufftPlanMany(&rfft_pad  , NRANK, dims_pad  , dims_pad     , 1, Nxy_pad     , dims_pad_half, 1, Nxy_pad_half, CUFFT_R2C, 1);
    cufftPlanMany(&rifft_pad , NRANK, dims_pad  , dims_pad_half, 1, Nxy_pad_half, dims_pad     , 1, Nxy_pad     , CUFFT_C2R, 1);
    cufftPlanMany(&ifft_probe, NRANK, dims_probe, dims_probe   , 1, Nxy_probe   , dims_probe   , 1, Nxy_probe   , CUFFT_C2C, 1);    
    float * d_slices;
    cudaMalloc( &d_slices,    num_slices_max*((long long) dimx)*dimy*sizeof(float) );
    cudaMemset(  d_slices, 0, num_slices_max*((long long) dimx)*dimy*sizeof(float) );
    float *d_slices_tmp;
    cudaMalloc( &d_slices_tmp, num_slices_max*Nxy*sizeof(float) );

    const mwSize size_slices[3]={dimx,dimy,num_slices_max};
    // plhs[6] = mxCreateNumericArray(3, size_slices, mxSINGLE_CLASS, mxREAL);
    // float * slices_before   = (float*)mxGetSingles(plhs[6]); 

    // plhs[7] = mxCreateNumericArray(3, size_slices, mxSINGLE_CLASS, mxREAL);
    // float * slices   = (float*)mxGetSingles(plhs[7]); 

    // plhs[8] = mxCreateNumericArray(3, size_slices, mxSINGLE_CLASS, mxREAL);
    // float * slices2   = (float*)mxGetSingles(plhs[8]); 

    // plhs[9] = mxCreateNumericArray(3, size_slices, mxSINGLE_CLASS, mxREAL);
    // float * slices3   = (float*)mxGetSingles(plhs[9]); 
    
    float2* d_float2[7];
    float*  d_float[4];
    //cudaMalloc( &temp, (9*Nxy_probe+3*Nxy_pad+6*Nxy_pad_half)*sizeof(float) );    
    cudaMalloc( &d_float2[0], Nxy_probe*sizeof(float2));
    cudaMalloc( &d_float2[1], Nxy_probe*sizeof(float2));
    cudaMalloc( &d_float2[2], Nxy_probe*sizeof(float2));
    cudaMalloc( &d_float2[3], Nxy_probe*sizeof(float2));
    cudaMalloc( &d_float[0] , Nxy_probe*sizeof(float));

    cudaMalloc( &d_float[1]  , Nxy_pad*sizeof(float));
    cudaMalloc( &d_float[2]  , Nxy_pad*sizeof(float));
    cudaMalloc( &d_float[3]  , Nxy_pad*sizeof(float));

    cudaMalloc( &d_float2[4] , Nxy_pad_half*sizeof(float2));
    cudaMalloc( &d_float2[5] , Nxy_pad_half*sizeof(float2));
    cudaMalloc( &d_float2[6] , Nxy_pad_half*sizeof(float2));  

    // iteration
    cout << "Starting reconstruction:" <<endl;
    for (long iter=0; iter<iterations; iter++){
        // set zero residual
        //setValue<<<blocksPerGridPrj, threadsPerBlock>>>(d_residual, 0.0, nPjsPoints);
        cudaMemset( d_residual, 0, nPjsPoints*sizeof(float) );

        // compute forward projection        
        for (int i = 0; i < Num_pjs; i++){            
            int k_min = pj_slices[2*i];
            int k_max = pj_slices[2*i+1];
            int num_slices = k_max - k_min +1 ;
            
            //const mwSize size_slices[3]={dimx,dimy,num_slices};
            //mxGPUArray * mex_slices = mxGPUCreateGPUArray(mwSize(3), size_slices, mxSINGLE_CLASS, mxREAL, MX_GPU_INITIALIZE_VALUES );            
            //float *d_slices = (float*) ( mxGPUGetData( mex_slices) );
            //RHS[0] = (mxArray*) mex_slices;
            //cout << "RHS[0]" <<endl;

            // radon transform for each d_slices
            RadonTF<<<blocksPerGridRec2, threadsPerBlock>>>(d_rec, d_Matrix + i*9, dimx,dimy,d_nc, dimx_ext,dimy_ext,d_nc_ext,
            o_ratio,d_xs2+i*o_ratio,d_ys2+i*o_ratio, d_slices, k_min, k_max, defocus_step );
            //d_residual + i*nrow_cols
            
            if(is_y_avg){                
                avg_y<<< dim3( (dimx+15)/16, (dimy+15)/16 ), dim3(16,16) >>>( d_slices, d_slices_tmp, dimx, dimy, num_slices);
                float *tmp   = d_slices;
                d_slices     = d_slices_tmp;
                d_slices_tmp = tmp;
            }

            //copy d_slices to RHS[0]:h_slices
            //const mwSize size_slices[3]={dimx,dimy,num_slices};
            //RHS[0] = mxCreateNumericArray(3, size_slices, mxSINGLE_CLASS, mxREAL);
            //float* h_slices = (float*)mxGetSingles(RHS[0]);
            
            //cudaMemcpy(h_slices, d_slices, dimx*dimy*num_slices*sizeof(float), cudaMemcpyDeviceToHost);

            // call Matlab function
            // const mwSize size_df_value[1] = {num_slices};
            // RHS[3] = mxCreateNumericArray(1, size_df_value, mxSINGLE_CLASS, mxREAL);
            // float * df_values   = (float*)mxGetSingles(RHS[3]);
            
            float * df_values = new float[num_slices];
            for(int k=k_min; k<=k_max; k++){
                df_values[k-k_min] = abs(k*defocus_step-defocus_param[i])*pixelsize * defocus_scale;
                //cout<< k <<"." << df_values[k-k_min]<<", ";
            }
            
            //mexCallMATLAB(1, LHS, 6, RHS, "FrequencyProbeFilter_forward"); 
            // plhs[1] = LHS[0];

            // get data from Matlab and copy to gpu
            // const float  * h_slices_blurry  = mxGetSingles(LHS[0]);
            // cudaMemcpy(d_slices, h_slices_blurry, num_slices*dimx*dimy*sizeof(float),  cudaMemcpyHostToDevice);


            /* to do: convolve d_slices with probe_filter
                for(int k=0, k<num_slices; k++){
                    d_slices[:,:,k] = conv2(d_slices[:,:,k], probe_filter) 
                }
            */

            int dims_slices[3] = { (int)dimx, (int)dimy, num_slices};

            //if(iter==iterations-1 && i==Num_pjs-1){      
                //cout << "num slices = " << num_slices << endl;          
                //cudaMemcpy(slices_before, d_slices, dimx*dimy*num_slices*sizeof(float), cudaMemcpyDeviceToHost);
                //const float  * h_slices_blurry  = mxGetSingles(LHS[0]);
                //memcpy(slices_before, h_slices_blurry, num_slices*dimx*dimy*sizeof(float) );
            //}


            //cout << "norm1 = " << sumGPUArray(d_slices,Nxy) <<endl; // (rfft_pad, rifft_pad, ifft_probe)
            //FrequencyProbeFilter_forward(d_slices, d_kR_rad2, d_ape, df_values, (float)lambda, hsize, dims_slices, dims_probe, rfft_pad, rifft_pad, ifft_probe);
            FrequencyProbeFilter(d_slices, d_kR_rad2, d_ape, df_values, (float)lambda, hsize, dims_slices, dims_probe, rfft_pad, rifft_pad, ifft_probe, d_float2, d_float);
            //cout << "norm2 = " << sumGPUArray(d_slices,Nxy) <<endl;
            
            //if(iter==iterations-1 && i==Num_pjs-1){                
                //cudaMemcpy(slices, d_slices, dimx*dimy*num_slices*sizeof(float), cudaMemcpyDeviceToHost);
            //}

            sum_slices<<<blocksPerGridSlice, threadsPerBlock>>>(d_residual + i*nrow_cols, d_slices, dimx,dimy , num_slices);
            
            cudaMemset(  d_slices, 0, num_slices_max*((long long) dimx)*dimy*sizeof(float) );
            //mxDestroyArray(RHS[3]);
            //mxDestroyArray(RHS[0]);
        }

        // copy calculated projections
        if(iter==iterations-1){
            cudaMemcpy( cal_proj, d_residual, Num_pjs*dimx*dimy*sizeof(float),  cudaMemcpyDeviceToHost);
        }

        // compute residual: = forward projection - measure projections
        computeResidual<<<blocksPerGridPrj, threadsPerBlock>>>(d_residual,d_projections,1.0f/o_ratio, nPjsPoints);

        // compute R1 factor
        if (iter%1 == 0){
            /*
            double res_norm, *res_norm_pt;
            cudaMalloc( (void**)&res_norm_pt, sizeof(double) );                
            cudaMemset(res_norm_pt, 0, 1*sizeof(double)); //setValue<<<1,1>>>(res_norm_pt, 0.0, 1);
            R1norm<<<blocksPerGridPrj, threadsPerBlock>>>(d_residual,    res_norm_pt, nPjsPoints);
            cudaMemcpy( &res_norm, res_norm_pt, sizeof(double), cudaMemcpyDeviceToHost ) ;
            cudaFree( res_norm_pt );  
            */
            // float res_norm, *sum_pt;            
            // cudaMalloc((void**)&sum_pt, sizeof(float)*gridSize);
            // sumCommMultiBlock<<<gridSize, blockSize>>>(d_residual, nPjsPoints, sum_pt);
            // sumCommMultiBlock<<<1, blockSize>>>(sum_pt, gridSize, sum_pt);            
            // //cudaDeviceSynchronize();
            // cudaMemcpy(&res_norm, sum_pt, sizeof(float), cudaMemcpyDeviceToHost);
            // cudaFree(sum_pt);

            float res_norm = sumGPUArray(d_residual, nPjsPoints);
            cout << iter+1 << ". R1 = " << res_norm/pj_norm << endl;
        }            

        // back projection
        for (long i = 0; i < Num_pjs; i ++){            
            int k_min = pj_slices[2*i];
            int k_max = pj_slices[2*i+1];
            int num_slices = k_max - k_min +1 ;
            
            // float * d_slices;
            // cudaMalloc( &d_slices,    num_slices*((long long) dimx)*dimy*sizeof(float) );

            //cudaMemset(  d_slices, 0, num_slices*((long long) dimx)*dimy*sizeof(float) );
            // to do: convolution d_slices with probe_filter instead of copy
            for(int k=0; k<num_slices; k++){
                cudaMemcpy( d_slices + k*dimx*dimy, d_residual + i*nrow_cols, (long long) dimx * dimy * sizeof(float), cudaMemcpyDeviceToDevice );            
            }
            

            //const mwSize size_slices[2]={dimx,dimy};
            //RHS[0] = mxCreateNumericArray(2, size_slices, mxSINGLE_CLASS, mxREAL);
            //float* h_slices = (float*)mxGetSingles(RHS[0]);
            //cudaMemcpy(h_slices, d_residual + i*nrow_cols, dimx*dimy*sizeof(float), cudaMemcpyDeviceToHost);

            // const mwSize size_df_value[1] = {num_slices};
            // RHS[3] = mxCreateNumericArray(1, size_df_value, mxSINGLE_CLASS, mxREAL);
            // float * df_values   = (float*)mxGetSingles(RHS[3]);
            
            float * df_values = new float[num_slices];
            for(int k=k_min; k<=k_max; k++){
                df_values[k-k_min] = abs(k*defocus_step-defocus_param[i])*pixelsize * defocus_scale;
            }

            //FrequencyProbeFilter_backward(A, kR_rad2, ape, df_series, lambda , hsize)
            //mexCallMATLAB(1, LHS, 6, RHS, "FrequencyProbeFilter_backward"); 
            
            // get data from Matlab and copy to gpu
            //const float  * h_slices_blurry  = mxGetSingles(LHS[0]);   
            //cudaMemcpy(d_slices, h_slices_blurry, ((long long) num_slices)*dimx*dimy*sizeof(float), cudaMemcpyHostToDevice);

            // call FrequencyProbeFilter_backward gpu
            int dims_slices[3] = { (int)dimx, (int)dimy, num_slices};

            // if(iter==iterations-1 && i==Num_pjs-1){                
            //     cudaMemcpy(slices2, d_slices       , dimx*dimy*num_slices*sizeof(float), cudaMemcpyDeviceToHost);
            // }

            //FrequencyProbeFilter_forward(d_slices, d_kR_rad2, d_ape, df_values, (float)lambda, hsize, dims_slices, dims_probe, rfft_pad, rifft_pad, ifft_probe);
            FrequencyProbeFilter(d_slices, d_kR_rad2, d_ape, df_values, (float)lambda, hsize, dims_slices, dims_probe, rfft_pad, rifft_pad, ifft_probe, d_float2, d_float);

            //if(iter==iterations-1 && i==Num_pjs-1){                
                //cudaMemcpy(slices3, d_slices       , dimx*dimy*num_slices*sizeof(float), cudaMemcpyDeviceToHost);
                //memcpy(slices3, h_slices_blurry, dimx*dimy*num_slices*sizeof(float) );
            //}

            RadonTpose_updateRec<<<blocksPerGridRec2, threadsPerBlock>>>(d_Matrix + i*9, dimx,dimy,d_nc, dimx_ext,dimy_ext,d_nc_ext, d_slices,
            o_ratio, d_xs2+i*o_ratio, d_ys2+i*o_ratio, d_rec, theta_dt[i], recPoints, k_min, k_max, defocus_step);
            

            cudaMemset(  d_slices, 0, num_slices_max*((long long) dimx)*dimy*sizeof(float) );
            //cudaFree(d_slices);
            //mxDestroyArray(RHS[3]);
            //mxDestroyArray(RHS[0]);
            //mxDestroyArray(LHS[0]);            
        }

        // apply constraint
        if(positivity && l2_reg>0){
            update_positivity_l2<<<blocksPerGridRec2, threadsPerBlock>>>(d_rec, dimx_ext, dimy_ext, scale_factor);
        }
        else if(positivity){
            update_positivity<<<blocksPerGridRec2, threadsPerBlock>>>(d_rec, dimx_ext, dimy_ext);
        }
        else if(l2_reg>0){
            update_l2<<<blocksPerGridRec2, threadsPerBlock>>>(d_rec, dimx_ext, dimy_ext, scale_factor);
        }

    }

    /* return result  */
    //plhs[0] = mxGPUCreateMxArrayOnCPU(rec);
    cufftDestroy(rfft_pad);
    cufftDestroy(ifft_probe);
    cufftDestroy(rifft_pad);
    
    cudaMemcpy( rec, d_rec, recPoints*sizeof(float), cudaMemcpyDeviceToHost);
    //plhs[1] = mxGPUCreateMxArrayOnGPU(residual);

    //mxGPUDestroyGPUArray(projections); 
    //mxGPUDestroyGPUArray(rec);        
    cudaFree(d_rec);
    cudaFree(d_residual);
    cudaFree(d_Matrix);
    cudaFree(d_projections);
    cudaFree(d_slices);
    cudaFree(d_slices_tmp);
    
    cudaFree(d_ape);
    cudaFree(d_kR_rad2);
    //mxGPUDestroyGPUArray(Matrix);
    //mxGPUDestroyGPUArray(nc);
    cudaFree(d_nc);
    //mxGPUDestroyGPUArray(xs2);
    //mxGPUDestroyGPUArray(ys2);
    //mxGPUDestroyGPUArray(residual);
    //mxGPUDestroyGPUArray(support); 
    cudaFree( shift_ptr );
    cudaFree( x_shift );
    cudaFree( y_shift );

    for(int d=0; d<7;d++)  cudaFree(d_float2[d]);
    for(int d=0; d<4;d++)  cudaFree(d_float[d]);    
    
    //cudaDeviceReset();      
    //mxGPUArray       *coord2 = (mxGPUArray*)d_Coord2;
    // plhs[1] = mxGPUCreateMxArrayOnCPU(coord2);
    
}


/*
// compute diffnorm of xyshift
float x_norm, x_norm2,  *x_pt, *x_pt2;
cudaMalloc( (void**)&x_pt, sizeof(float) );
cudaMalloc( (void**)&x_pt2, sizeof(float) );
cudaMemset(x_pt, 0, sizeof(float));
cudaMemset(x_pt2, 0, sizeof(float));
R1norm<<<2, threadsPerBlock>>>(d_xs,  x_pt, Num_pjs*4);
R1norm<<<2, threadsPerBlock>>>(d_xs2,  x_pt2, Num_pjs*4);
cudaMemcpy( &x_norm ,  x_pt, sizeof(float), cudaMemcpyDeviceToHost ) ;
cudaMemcpy( &x_norm2 ,  x_pt2, sizeof(float), cudaMemcpyDeviceToHost ) ;
cudaFree( x_pt);
cout << "x_norm = " << x_norm << ", x_norm2 = " << x_norm2 <<endl;
*/




















//
//  con_level.cpp
//
//
//  Created by aylab on 7/18/16.
//
//
#ifndef CONCENTRATION_LEVEL_HPP
#define CONCENTRATION_LEVEL_HPP
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <iostream>
#include "macros.hpp"

#ifdef __CUDACC__
#define CPUGPU_FUNC __host__ __device__
#else
#define CPUGPU_FUNC
#endif

template<class T>
class concentration_level{
public:
    
    class cell{
    public:
        CPUGPU_FUNC cell(T *row): _array(row) {}
        CPUGPU_FUNC T& operator[](int k){
			return _array[k];
		}
        CPUGPU_FUNC const T& operator[](int k) const {
			return _array[k];
		}
        T *_array;
    };


    class timespan{
    public:
        CPUGPU_FUNC timespan(T *plane,int width): _array(plane), _width(width) {};
        CPUGPU_FUNC cell operator[](int j) {
            cell temp(_array+_width*j);
        	return temp;
        }

        CPUGPU_FUNC const cell operator[](int j) const{
            cell temp(_array+_width*j);
        	return temp;
        }
        T *_array;
        int _width;
    };
    concentration_level(int height =0, int length =0, int width =0)
		:_height(height),_length(length),_width(width),_cuda(false){
        allocate_array();
    }
 
#if 0   
    concentration_level(const concentration_level<T>& other)
		:_height(other._height),_length(other._length),_width(other._width),_cuda(other._cuda){
        allocate_array();
        for (int he=0; he< _height; he++){
            for (int le=0; le< _length; le++){
                for (int wi=0; wi< _width; wi++){
                    _array[he*_length*_width+le*_width+wi]=other._array[he*_length*_width+le*_width+wi];
					_darray[he*_length*_width+le*_width+wi]=other._darray[he*_length*_width+le*_width+wi];
                }
            }
        }
    }
#endif
    
    void initialize(int height, int length, int width){
        dealloc_array();
        _width=width;
        _height=height;
        _length=length;
		_cuda=false;
        allocate_array();
        reset();
    }
    
    void reset(){
        for (int i = 0; i < _height; i++) {
            for (int j = 0; j < _length; j++) {
                for (int k = 0; k < _width; k++) {
                    _array[i*_length*_width+j*_width+k] = 0; // Initialize every concentration level at every time step for every cell to 0
					
                }
            }
        }
    }
    
	/*
	concentration_level& operator=(const concentration_level& other){
        dealloc_array();
        _height=other._height;
        _length=other._length;
        _width=other._width;
        allocate_array();
        for (int he; he< _height; he++){
            for (int le; le< _length; le++){
                for (int wi; wi< _width; wi++){
                    _array[he*_length*_width+le*_width+wi]=other._array[he*_length*_width+le*_width+wi];
                }
            }
        }
        return *this;
    }
    */
	CPUGPU_FUNC timespan operator[](int i){
		if (_cuda){
			timespan temp(_darray+_length*_width*i, _width);
       		return temp;
		}		
		else{
			timespan temp(_array+_length*_width*i, _width);
        	return temp;
		}		
	}
	
    CPUGPU_FUNC const timespan operator[](int i) const{
		if (_cuda){
			timespan temp(_darray+_length*_width*i, _width);
       		return temp;
		}		
		else{
			timespan temp(_array+_length*_width*i, _width);
        	return temp;
		}		
	}
/*    
    ~concentration_level() {
        dealloc_array();
    }
 */   
    int height() const {return _height;}
    int length() const {return _length;}
    int width() const {return _width;}
    CPUGPU_FUNC bool getStatus() { return _cuda; }

	void toGPU(){
		if(_cuda){
			return;
		}	
		//std::cout << "Copying to GPU\n";
		int size = _height*_length*_width* sizeof(T);
		CUDA_ERRCHK(cudaMalloc((void**)&_darray, size));
		CUDA_ERRCHK(cudaMemcpy(_darray,_array,size,cudaMemcpyHostToDevice));
		_cuda=true;
	}

	void toCPU(){
		if(!_cuda){
			return;
		}
		//std::cout << "Copying to CPU\n";
		int size = _height*_length*_width* sizeof(T);
		CUDA_ERRCHK(cudaMemcpy(_array, _darray, size, cudaMemcpyDeviceToHost));
		CUDA_ERRCHK(cudaFree(_darray));
		_darray = NULL;
		_cuda= false;
	}
protected:
    void dealloc_array(){
        if (_array){
            delete[] _array;
        }
        _array= NULL;
    }
	
    void allocate_array(){
        if (_width*_length*_height >0){
            _array= new T[_height*_length*_width];
			if (_array == NULL){std::cout<<"ERROR"<<std::endl; exit(EXIT_MEMORY_ERROR);}
        }
        else{
            _array= NULL;
		}
    }
    
    
    int   _height,_length, _width;
	bool _cuda;
    T *_array;
	T *_darray;
};
#endif

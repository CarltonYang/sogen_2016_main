//
//  con_level.cpp
//
//
//  Created by aylab on 7/18/16.
//
//
#ifndef ARRAY2D_HPP
#define ARRAY2D_HPP
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
class array2D{
public:
    
    class cell{
    public:
        CPUGPU_FUNC cell(T *row, bool cuda): _array(row), _cuda(cuda) {}
        CPUGPU_FUNC T& operator[](int k){
			if (_cuda){return _darray[k];}
			else{return _array[k];}
		}
        CPUGPU_FUNC const T& operator[](int k) const {
			if (_cuda){return _darray[k];}
			else{return _array[k];}
		}
        T *_array;
		T *_darray;
		bool _cuda;
    };


   
    array2D(int height =0,  int width =0)
		:_height(height),_width(width),_cuda(false){
        allocate_array();
    }
    
    array2D(const array2D<T>& other)
		:_height(other._height),_width(other._width),_cuda(other._cuda){
        allocate_array();
        for (int he=0; he< _height; he++){
                for (int wi=0; wi< _width; wi++){
                    _array[he*_width+wi]=other._array[he*_width+wi];
					_darray[he*_width+wi]=other._darray[he*_width+wi];
                }
            }
        }
    }
    
    void initialize(int height, int width){
        dealloc_array();
        _width=width;
        _height=height;
		_cuda=false;
        allocate_array();
        reset();
    }
    
    void reset(){
        for (int i = 0; i < _height; i++) {
                for (int k = 0; k < _width; k++) {
                    _array[i*_width+j*_width+k] = 0; // Initialize every concentration level at every time step for every cell to 0
					
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
	CPUGPU_FUNC cell operator[](int i){
		if (_cuda){
			cell temp(_darray+_width*i,  _cuda);
       		return temp;
		}		
		else{
			cell temp(_array+_width*i,  _cuda);
        	return temp;
		}		
	}
	
    CPUGPU_FUNC const timespan operator[](int i) const{
		if (_cuda){
			cell temp(_darray+_width*i, _cuda);
       		return temp;
		}		
		else{
			cell temp(_array+_width*i, _cuda);
        	return temp;
		}		
	}
    
    ~array2D() {
        dealloc_array();
    }
    
    int height() const {return _height;}
    int width() const {return _width;}
    CPUGPU_FUNC bool getStatus() { return _cuda; }

	void toGPU(){
		if(!_cuda){
			return;
		}	
		int size = _height*_width* sizeof(T);
		cudaMalloc((void**)&_darray, size);
		cudaMemcpy(_darray,_array,size,cudaMemcpyHostToDevice);
		_cuda=true;
	}

	void toCPU(){
		if(!_cuda){
			return;
		}
		int size = _height*_width* sizeof(T);
		cudaMemcpy(_array, _darray, size, cudaMemcpyDeviceToHost);
		cudaFree(_darray);
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
        if (_width*_height >0){
            _array= new T[_height*_width];
			if (_array == NULL){std::cout<<"ERROR"<<std::endl; exit(EXIT_MEMORY_ERROR);}
        }
        else{
            _array= NULL;
		}
    }
    
    
    int   _height, _width;
	bool _cuda;
    T *_array;
	T *_darray;
};
#endif
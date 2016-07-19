//
//  con_level.cpp
//
//
//  Created by aylab on 7/18/16.
//
//
#ifndef CONCENTRATION_LEVEL_HPP
#define CONCENTRATION_LEVEL_HPP

template<class T>
class concentration_level{
public:
    
    class cell{
    public:
        cell(T *row): _array(row) {}
        T& operator[](int k){return _array[k];}
        const T& operator[](int k) const {return _array[k];}
        T *_array;
    };
    class timespan{
    public:
        timespan(T *plane,int width): _array(plane), _width(width){};
        cell operator[](int j) {
            cell temp(_array+_width*j);
            return temp;
        }
        const cell operator[](int j) const{
            cell temp(_array+_width*j);
            return temp;
        }
        T *_array;
        int _width;
    };
    concentration_level(int height =0, int length =0, int width =0):_height(height),_length(length),_width(width){
        allocate_array();
    }
    
    concentration_level(const concentration_level<T>& other):_height(height),_length(length),_width(width){
        allocate_array();
        for (int he; he< _height; he++){
            for (int le; le< _length; le++){
                for (int wi; wi< _width; wi++){
                    _array[he*_length*_width+le*_width+wi]=other._array[he*_length*_width+le*_width+wi];
                }
            }
        }
    }
    
    void initialize(int height, int length, int width){
        dealloc_array();
        _width=width;
        _height=height;
        _length=length;
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
    
    timespan operator[](int i){
        timespan temp(_array+_length*_width*i, _width);
        return temp;
    }
    
    const timespan operator[](int i) const{
        timespan temp(_array+_length*_width*i, _width);
        return temp;
    }
    
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
    
    ~concentration_level() {
        dealloc_array();
    }
    
    int height() const {return _height;}
    int length() const {return _length;}
    int width() const {return _width;}
    
protected:
    void dealloc_array(){
        if (_array){
            delete[] _array;
        }
        _array= NULL;
    }
    
    void allocate_array(){
        if (_width*_length*_height >0){
            _array= new T[_width*_length*_height];
        }
        else
            _array= NULL;
    }
    
    
    int   _height,_length, _width;
    T *_array;
};
#endif
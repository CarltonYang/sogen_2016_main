//
//  con_level.cpp
//  
//
//  Created by aylab on 7/18/16.
//
//

#include "con_level.hpp"

template<class T>
class con_levels{
public:
    class timespan{
    public:
        timespan(int length, int width): _length(length), _width(width);
        cell& operator[](int j) {
            cell temp(_array+_width*j);
            return temp;
        }
        const cell& operator[](int j) const{
            cell temp(_array+_width*j);
            returm temp;
        }
    };
    class cell{
    public:
        cell(T *cell): _array(cell) {}
        T& operator[](int k){return _array[j];}
        const T& operator[](int k) const {return _array[j];}
        
        T *_array;
    };
    
    con_levels(int height =0, int length =0, int width =0):_height(height),_length(length),_width(width){
        allocate_array();
    }
    
    con_levels(const con_levels<T>& other):_height(height),_length(length),_width(width){
        allocate_array();
        for (int he; he< _height; he++){
            for (int le; le< _length; le++){
                for (int wi; wi< _width; wi++){
                    _array[he*_length*_width+le*_width+wi]=other._array[he*_length*_width+le*_width+wi];
                }
            }
        }
    }
    timespan operator[](int i){
        timespan temp(_array+_length*_width*i);
        return temp;
    }
    
    const timespan operator[](int i) const{
        timespan temp(_array+_length*_width*i);
        return temp;
    }
    
    con_levels& operator=(const con_levels& other){
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
    
    ~con_levels() {
        dealloc_array();
    }
};

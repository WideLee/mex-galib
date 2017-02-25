//
// Created by amax on 17-2-15.
//

#ifndef GATEST_MATRIX_H
#define GATEST_MATRIX_H
//
//  iBS_Matrix.h
//
//
//  Created by nash on 11/29/15.
//  Copyright 2015 iBean Software.
//  All rights reserved.
//

const int Matrix_MAJOR_VERSION = 1;
const int Matrix_MINOR_VERSION = 0;

#include <iostream>
#include <vector>
namespace iBS
{
    struct Matrix
    {
        std::vector<std::vector<double> > data;

        Matrix(int size) {
            data.resize(size);
            for (int i = 0; i < size; ++i) {
                data[i].resize(size);
            }
        }

        Matrix& operator =(Matrix& o)
        {
            data.resize(o.data.size());
            for(int i=0;i<data.size();i++)
                data[i].resize(o.data[i].size());
            for(int i=0;i<data.size();i++)
                for(int j=0;j<data[i].size();j++)
                {
                    data[i][j] = o.data[i][j];
                }
            return *this;
        }

        Matrix& operator +(Matrix& o)
        {
            for(int i=0;i<data.size();i++)
                for(int j=0;j<data[i].size();j++)
                {
                    data[i][j] = data[i][j] + o.data[i][j];
                }
            return *this;
        }
        Matrix& operator -(Matrix& o)
        {
            for(int i=0;i<data.size();i++)
                for(int j=0;j<data[i].size();j++)
                {
                    data[i][j] = data[i][j] - o.data[i][j];
                }
            return *this;
        }
        Matrix& operator *(Matrix& o)
        {
            if(data[0].size() != o.data.size()) return *this;

            Matrix tm(1);
            tm.data.resize(data.size());
            for(int i=0;i<tm.data.size();i++)
                tm.data[i].resize(o.data[0].size());

            for(int i=0;i<tm.data.size();i++)
                for(int j=0;j<tm.data[i].size();j++)
                {
                    tm.data[i][j] = 0;
                    for (int c=0; c<data[i].size(); c++)
                    {
                        tm.data[i][j] += data[i][c] * o.data[c][j];
                    }

                }
            *this = tm;
            return *this;
        }
        Matrix& operator ^(int power)
        {
            Matrix  tM2(1);
            tM2 = *this;

            //   not <= below \/ because first time counts as 2
            for(int i=1; i<power; ++i)
                *this = (*this) * (tM2);

            return *this;
        }

        void print()
        {
            for(int i=0;i<data.size();i++)
            {
                for(int j=0;j<data[i].size();j++)
                {
                    std::cout << data[i][j] << ' ';
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    };

}; // end of namespace iBS

#endif //GATEST_MATRIX_H

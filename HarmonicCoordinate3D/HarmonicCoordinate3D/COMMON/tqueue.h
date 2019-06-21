//--------------------------------------------------------------------------
// This file is released under 3-clause modified bsd license.
//
// Copyright (c) 2016, Takashi Ijiri
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//* Redistributions of source code must retain the above copyright notice, 
//  this list of conditions and the following disclaimer.
//* Redistributions in binary form must reproduce the above copyright notice, 
//  this list of conditions and the following disclaimer in the documentation 
//  and/or other materials provided with the distribution.
//* Neither the names of the Ritsumeikan University nor the names of its contributors 
//  may be used to endorse or promote products derived from this software 
//  without specific prior written permission.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
//DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
//(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
//LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ---------------------------------------------------------------------------



/* -----------------------------------------------------------------
 * queue to replace stl::queue and stl::deque
 * Simple and fast since this class perform no error check
 * "head" indicates the head element
 * "tail" indicates the next of the last element (常に一つは空白がある)
 *
 * CG Gems JP 2013/2014用にOpenSourceとして開発を開始
 * 2015/9/7 に再メンテナンス
-------------------------------------------------------------------*/
#pragma once

#include <iostream>

#pragma unmanaged

#define T_QUEUE_INIT_SIZE 40
#define T_QUEUE_ADD_SIZE  1000

template<class T>
class TQueue
{
  const int m_INCREASE_SIZE;
  int m_size, m_tail, m_head;
  T  *m_data;
public:
  ~TQueue() { delete[] m_data; }

  //init_size     : キューが最初に用意しておくメモリサイズ (要素数)
  //increase_size : 要素数が用意したメモリの上限に達したときに再確保するメモリサイズ (要素数)
  TQueue(const int init_size = T_QUEUE_INIT_SIZE, const int increase_size = T_QUEUE_ADD_SIZE) : m_INCREASE_SIZE( std::max(20, increase_size))
  {
    std::cout << "TQueue standard constractor\n";
    m_size = std::max(10, init_size);
    m_data = new T[m_size];
    m_tail = m_head = 0;
  }

  TQueue(const TQueue &src) : m_INCREASE_SIZE(src.m_INCREASE_SIZE) {
    std::cout << "TQueue copy constractor\n";
    m_size = src.m_size;
    m_tail = src.m_tail;
    m_head = src.m_head;
    m_data = new T[m_size];
    //memcpy( m_data, src.m_data, sizeof( T ) * m_size ); ← これでは浅いコピーになりエラー
    if (m_head < m_tail) for (int i = m_head; i < m_tail; ++i) m_data[i] = src.m_data[i];
    else {
      for (int i = m_head; i < m_size; ++i) m_data[i] = src.m_data[i];
      for (int i = 0; i < m_tail; ++i) m_data[i] = src.m_data[i];
    }
  }

  TQueue& operator=(const TQueue& src) {
    std::cout << "TQeue operator = \n";
    delete[] m_data;
    m_size = src.m_size;
    m_tail = src.m_tail;
    m_head = src.m_head;
    m_data = new T[m_size];
    //memcpy( m_data, src.m_data, sizeof( T ) * m_size ); ← これでは浅いコピー
    if (m_head < m_tail) for (int i = m_head; i < m_tail; ++i) m_data[i] = src.m_data[i];
    else {
      for (int i = m_head; i < m_size; ++i) m_data[i] = src.m_data[i];
      for (int i = 0; i < m_tail; ++i) m_data[i] = src.m_data[i];
    }
    return *this;
  }

  inline void swap(TQueue &trgt) {
    int t; T* p;
    t = m_size; m_size = trgt.m_size; trgt.m_size = t;
    t = m_tail; m_tail = trgt.m_tail; trgt.m_tail = t;
    t = m_head; m_head = trgt.m_head; trgt.m_head = t;
    p = m_data; m_data = trgt.m_data; trgt.m_data = p;
  }

  inline bool empty() { return m_tail == m_head; }
  inline bool hasElement() { return m_tail != m_head; }
  inline void clear() { m_tail = m_head = 0; }
  inline T&   front() { return m_data[m_head]; }
  inline T&   back() { return m_data[m_tail - 1]; }
  inline void pop_front() { m_head++; if (m_head == m_size) m_head = 0; }
  inline void pop_back() { m_tail--; if (m_tail < 0) m_tail = m_size - 1; }

  inline int  size()const { //要素数を返す 
    return (m_tail >= m_head) ? m_tail - m_head : m_tail + m_size - m_head;
  }

  //i番目の要素を返す
  inline   T& operator[](const int &i) { int idx = i + m_head; if (idx >= m_size) idx -= m_size; return m_data[idx]; }
  inline   T  operator[](const int &i) const { int idx = i + m_head; if (idx >= m_size) idx -= m_size; return m_data[idx]; }

  inline void push_back(const T &n)
  {
    m_data[m_tail] = n;
    m_tail++;
    if (m_tail == m_size) m_tail = 0;
    if (m_tail == m_head) increase_size();//サイズ拡張
  }
  inline void push_front(const T &n)
  {
    --m_head;
    if (m_head < 0) m_head = m_size - 1;
    m_data[m_head] = n;
    if (m_tail == m_head) increase_size();//サイズ拡張
  }

private:
  inline void increase_size()
  {
    std::cout << "*";
    T* newData = new T[m_size + m_INCREASE_SIZE];

    //head -- m_size-1 をコピー (memcpyは浅いコピーだから使っちゃダメ)
    int newI = 0;
    for (int i = m_head; i < m_size; ++i, ++newI) newData[newI] = m_data[i];
    for (int i = 0; i < m_tail; ++i, ++newI) newData[newI] = m_data[i];

    delete[] m_data;
    m_data = newData;

    m_head = 0;
    m_tail = m_size;
    m_size = m_size + m_INCREASE_SIZE;
  }

public:
  //debug用
  void vis()
  {
    std::cout << "\n";
    std::cout << m_size << " " << m_head << " " << m_tail << "\n";
    for (int i = 0; i < m_size; ++i) {
      if (i == m_head) std::cout << "h";
      if (i == m_tail) std::cout << "t";
      std::cout << m_data[i] << " ";
    }
    std::cout << "\n";
  }
  static void test() {

    TQueue<int> Q(5, 5);
    for (int i = 3; i > -10; --i) { Q.push_front(i); Q.vis(); }
    for (int i = 3; i < 15; ++i) { Q.push_back(i); Q.vis(); }
    Q.pop_front();
    Q.pop_front();
    Q.pop_front();
    Q.pop_front();
  }
};


#pragma managed


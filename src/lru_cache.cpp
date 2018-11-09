/*
 * 
 *
 * File:   lrucache.hpp
 * Author: Alexander Ponomarev
 *
 * Created on June 20, 2013, 5:09 PM
 * 
 * Copyright (c) 2014, lamerman
 * All rights reserved.
 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
* Neither the name of lamerman nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.
*THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#ifndef _LRUCACHE_HPP_INCLUDED_
#define	_LRUCACHE_HPP_INCLUDED_

#include <unordered_map>
#include <list>
#include <cstddef>
#include <stdexcept>
#include <RcppArmadillo.h>

namespace cache {

//Cache stores a vector of type value_t for each input.
template<typename key_t, typename value_t>
class lru_cache_vec {
public:
  typedef typename std::pair<key_t, size_t> key_index_pair_t;
  typedef typename arma::Mat<value_t> cache_t;
  typedef typename arma::Col<value_t> value_vec_t;
  typedef typename std::list<key_index_pair_t>::iterator list_iterator_t;
  
  lru_cache_vec(size_t max_size, size_t value_length) :
    _max_size(max_size), _value_length(value_length), 
    _cache_mem(value_length, max_size) {
    }
  
  void put(const key_t& key, const value_vec_t& val_vec) {
    
    if(_max_size < 1) return; //Nothing to do
    
    auto it = _cache_items_map.find(key);
    if (it != _cache_items_map.end()) {
      //Key already in cache. 
      //Just update the corresponding memory location and realign the queue
      _cache_items_list.splice(_cache_items_list.begin(), _cache_items_list, it->second);
      _cache_mem.col(it->second->second) = val_vec;
      return;
    }
    
    //Now _max_size >= 1 and item is not already in the cache.
    
    //The index to store the val_vec at
    size_t next_index = _cache_items_map.size();
    
    if (next_index >= _max_size) {
      //cache full. First create some space
      auto last = _cache_items_list.end();
      last--;
      
      //The freed-up index
      next_index = last->second;
      
      //Delete the old_item from map and list
      _cache_items_map.erase(last->first);
      _cache_items_list.pop_back();
    }
    
    //Add the new key to map and list.
    _cache_items_list.push_front(key_index_pair_t(key, next_index));
    _cache_items_map[key] = _cache_items_list.begin();
    
    //Finally update the memory at that index
    _cache_mem.col(next_index) = val_vec;
  }
  
  value_vec_t get(const key_t& key) {
    auto it = _cache_items_map.find(key);
    if (it == _cache_items_map.end()) {
      throw std::range_error("There is no such key in cache");
    } else {
      _cache_items_list.splice(_cache_items_list.begin(), _cache_items_list, it->second);
      return _cache_mem.col(it->second->second);
    }
  }
  
  bool exists(const key_t& key) const {
    return _cache_items_map.find(key) != _cache_items_map.end();
  }
  
  size_t size() const {
    return _cache_items_map.size();
  }
  
private:
  std::list<key_index_pair_t> _cache_items_list;
  std::unordered_map<key_t, list_iterator_t> _cache_items_map;
  size_t _max_size, _value_length;
  cache_t _cache_mem;
};

} // namespace cache

#endif	/* _LRUCACHE_HPP_INCLUDED_ */
import pickle
from element import element
from typing import List

def read_element(filename:str) -> List[element]:
  '''
  此函数从二进制文件中读取序列化的单元类列表
  filename: 二进制文件的名称
  返回值: 无
  '''
  with open(file=filename, mode='rb') as f:
    my_restored_list = pickle.load(f)
  return my_restored_list

element_list = read_element("./tmp.pickle")
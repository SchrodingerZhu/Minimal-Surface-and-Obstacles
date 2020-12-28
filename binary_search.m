%%
%% Binary Linear Searcher (used in PW condition search)
%% Created by Schrodinger ZHU <i@zhuyi.fan> for DDA/CIE6010 Project
%%
%% USAGE:
%% 1.     l, lowerbound
%% 2.     u, upperbound
%% 3. check, checking function return   
%% 4.   eps, error bound
%%

function [m] = binary_search(l, u, check, eps)
  m = (u + l) / 2;
  res = check(m);
  while u - l > eps && res ~= 0
      if res > 0
          u = m;
      else
          l = m;
      end
      m = (u + l) / 2;
      res = check(m);
  end
end

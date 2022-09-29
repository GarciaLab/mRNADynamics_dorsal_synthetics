function [idx, val] = nearestIndex(v, q)

   [~, idx] = min(abs(v-q));
   val = v(idx);
    
end
function [result] = B_hat(A,B,cur_state,final_state)

del = final_state - cur_state;

for i = del : -1 : 1
    if i == del
        result(:,:,i) = B(:,:,final_state-1);
    else
        result(:,:,i) = mA(A,cur_state+i,final_state-1) * B(:,:,cur_state + i -1);
    end    
end

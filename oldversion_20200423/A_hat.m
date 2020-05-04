function [result] = A_hat(A,cur_state,final_state)

[result] = mA(A,cur_state, final_state-1);

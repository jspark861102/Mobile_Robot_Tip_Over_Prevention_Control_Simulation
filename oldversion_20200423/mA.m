function [multiple_A] = mA(A,cur_state, final_state)

multiple_A = eye(size(A(:,:,1)));
for i = cur_state : final_state
    multiple_A = A(:,:,i)*multiple_A;
end
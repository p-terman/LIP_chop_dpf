% This one translate the hit vector in xlm to readable form.
function HitVector = XLMHitVectorReadout(hit_vector)
    % Ensure data type
    hit_vector_unsigned = uint8(hit_vector);
    BinaryS1HitVector = dec2bin(hit_vector_unsigned);%Check the number from base 10 to base 2. The output is a string
    [r c] = size(BinaryS1HitVector);%Check the original size of BinaryS1HitVector, hence, r should be twice of total events loaded, since dec2bin returns linear index form, and c would be for 8 channel
    HitFlagTemp = zeros(r,8);%Initialize matrix for keeping hit flag
    for i=1:r
        for j=1:c%Read the hit flag, change them to number.
            HitFlagTemp(i,j)=str2double(BinaryS1HitVector(i,c-j+1));
        end
    end
    HitVector = horzcat( HitFlagTemp(1:(r/2),:) , HitFlagTemp(((r/2)+1):end,:) ) ;
end
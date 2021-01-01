function C = nstruct2cell(S, BRANCH)
%C = NSTRUCT2CELL(S)
%
% Recursive function that converts a nested struct S with a total of n sub-fields into a nx2 cell array C.
% The first column of C lists the full names of the sub-fields and the second column contains the respective content.

% Mathias Benedek 201-01-04


if nargin == 1  %First level of struct
    BRANCH = inputname(1);
end
C = {};

if ~isstruct(S)     % End of struct-branch, no further fields: read content
    C = {BRANCH, S};

else
    n = numel(S);

    if n == 1   % (non-array) struct: parse fields
        fn = fieldnames(S);
        for ii = 1:length(fn)
            C = [C; nstruct2cell(S.(fn{ii}), [BRANCH,'.',fn{ii}])];
        end

    else        % struct-array: parse array elements
        for jj = 1:n
            C = [C; nstruct2cell(S(jj), [BRANCH,'(',num2str(jj),')'])];
        end

    end

end
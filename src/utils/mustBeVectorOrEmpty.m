function mustBeVectorOrEmpty(x)
%MUSTBEVECTOROREMPTY Check if x is empty or a vector

if ~isvector(x) && ~isempty(x)
    error("Input must be vector or empty.");
end

end


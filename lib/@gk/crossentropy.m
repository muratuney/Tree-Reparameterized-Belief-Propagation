function d = crossentropy( a, b)

if ~isa( a, 'gk' )
    error('The first argument should be a gk object');
else
    if length(a)>1
        error('The first argument should be a single object');
    end
end

if isa( b, 'gk')
    c = b;
else
    error('The second argument should be an array of gk objects.')
end

d = zeros(size(b));
for i=1:length(b(:))
    d(i) = kld(a, b(i)) + entropy( a );
end

function [ret, ret2] = eval_seg_ggs_divide( data, seg_point, lambda )
[T, col] = size( data );
reg_mat = eye( col ) * lambda;
chunk_size = length(seg_point)+1;
ret2 = zeros( chunk_size, 1 );
for i=1:chunk_size
    if i==1
        tdata = data( 1:(seg_point(i)-1), : );
    elseif i == chunk_size
        tdata = data( seg_point(i-1):T, : );
    else
        tdata = data( seg_point(i-1):(seg_point(i)-1), : );
    end
    t_size = size(tdata, 1);
    if t_size == 1
        cov_mat = reg_mat;
    else
        cov_mat = cov( tdata ) + ( reg_mat / t_size );
    end
    if sum(sum(isnan( cov_mat ))) ~= 0 || sum(sum(isinf( cov_mat ))) ~= 0
        fprintf('tdata size: %d %d\n', size(tdata, 1), size(tdata, 2));
        fprintf('i: %d seg_point(i): %d\n', i, seg_point(i));
    end
    ret2(i, 1) = -0.5 * t_size * log( det( cov_mat ) ) - lambda * trace( pinv(cov_mat) );
end
ret = sum(ret2);
end

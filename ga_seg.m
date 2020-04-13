% reorder codes
function [ret, total_res] = ga_seg( data, pop_size, gen_size, gen_num, pc, pm, max_bp, lambda )
%total_res = zeros( max_bp + 1, 1 );
%% RANDOM INITIALIZE
T = size( data, 1 );
pop = zeros( pop_size, max_bp );
unique_size = 0;
while unique_size ~= pop_size
    for i=(unique_size+1):pop_size
        t1 = randperm( T-1 ) + 1;
        pop( i, : ) = sort(t1(1:max_bp));
    end
    pop = unique( pop, 'rows' );
    unique_size = size( pop, 1 );
end

%% EVALUATE
eval_val = zeros( pop_size, 1 );
for i=1:pop_size
	eval_val(i) = eval_seg_ggs( data, pop(i), lambda );
end
[eval_val, idx] = sort( eval_val, 'descend' );
fprintf('best: %.2f\n', eval_val(1, 1));
total_res( 1, 1 ) = eval_val(1, 1);
pop = pop( idx, : );

%% ITERATION
for gg=1:gen_num
    %% CROSSOVER
    num_cross = round(gen_size * pc);
    num_mu = round(gen_size * pm);
    
    offs = zeros( num_cross + num_mu, max_bp );
    for i=1:num_cross
        t1 = randperm( round(pop_size/10) );
        chr1 = pop( t1(1), : );
        chr2 = pop( t1(2), : );
        offs(i, :) = crossover( chr1, chr2, max_bp );
    end

    %% MUTATION
    for i=(num_cross+1):(num_cross+num_mu)
        offs(i, :) = mutation( pop( randi(round(pop_size/10)), : ), T, max_bp );
    end

    %% FILTER
    offs = filter_pop( pop, offs );
    new_eval = zeros( size(offs, 1), 1 );
    for i=1:size(offs, 1)
        new_eval(i) = eval_seg_ggs( data, offs(i, :), lambda );
    end
    
    %% SELECTION
    merge_pop = [pop; offs];
    [merge_eval, sort_idx] = sort( [eval_val; new_eval], 'descend' );
    pop = merge_pop( sort_idx(1:pop_size), : );
    eval_val = merge_eval(1:pop_size);
    total_res( gg+1, 1 ) = eval_val(1, 1);
    
    %% PRINT
    if mod( gg, 10 ) == 0
        fprintf('best: %.2f\n', eval_val(1, 1));
    end
end

ret = pop;

end	% function seg_ga end

function ret = crossover( chr1, chr2, max_bp )
% ret vitally contain redundant point
idx1 = ismember( chr1, chr2 );
idx2 = ismember( chr2, chr1 );
red_point = chr1( idx1 );
if isempty( red_point )
    t1 = randperm( max_bp-1 );
    ret = sort( [chr1( 1: t1 ), chr2( (t1+1):end )] );
elseif length( red_point ) == max_bp
    ret = chr1;
else
    chr1( idx1 ) = [];
    chr2( idx2 ) = [];
    uni_val = unique( [chr1, chr2] );
    num_remain = max_bp - length(red_point);
    t1 = randperm( length(uni_val) );
    ret = sort( [uni_val( t1(1:num_remain) ), red_point] );
end
end

function ret = mutation( chr, T, max_bp )
col = length( chr );
t1 = randi( col );
chr( t1 ) = randi( T-1 ) + 1;
while length(unique(chr)) ~= max_bp
    t1 = randi( col );
    chr( t1 ) = randi( T-1 ) + 1;
end
ret = sort(chr);
end

function ret = filter_pop( pop, new_gen )
pop_size = size( pop, 1 );
gen_size = size( new_gen, 1 );
idx = zeros( 0, 1 );
for i=1:gen_size
	check_same = 0;
	for j=1:pop_size
		if isequal( new_gen(i, :), pop(j, :) )
			check_same = 1;
			break;
		end
	end
	if check_same == 0
		idx( end+1, 1 ) = i;
	end
end
ret = new_gen( idx, : );
end

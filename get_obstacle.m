
function [X, Y, Z] = get_obstacle(x, y, div, kind, opts)
    if strcmp(kind, 'rect')
        [X, Y] = meshgrid(x:div:x+opts.length, y:div:y+opts.width);
        Z = opts.height * ones(size(X));
        X = reshape(X, [], 1);
        Y = reshape(Y, [], 1);
        Z = reshape(Z, [], 1);
    elseif strcmp(kind, 'semisphere')
        R = opts.radius;
        X = [];
        Y = [];
        for i = x-R:div:x+R
            for j = y-R:div:y+R
                if (i - x)^2  + (j - y)^2 <= R^2
                    X(end+1,1) = i;
                    Y(end+1,1) = j;
                end
            end
        end
        Z = sqrt(R^2 - (X-x).^2 - (Y-y).^2) / div;
    elseif strcmp(kind, 'cone')
        R = opts.radius;
        H = opts.height;
        X = [];
        Y = [];
        Z = [];
        for i = x-R:div:x+R
            for j = y-R:div:y+R
                if (i - x)^2  + (j - y)^2 <= R^2
                    X(end+1,1) = i;
                    Y(end+1,1) = j;
                end
            end
        end
        Z = H .* (R - sqrt((X-x).^2 + (Y-y).^2)) ./ R;
    end
    X = round(X/div);
    Y = round(Y/div);
end
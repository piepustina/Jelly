classdef MyHandleClass < handle & my.super.Class
    % a handle class
    %
    % :param x: a variable

    %% some comments
    properties
        x % a property

        % Multiple lines before a
        % property can also be used
        y
    end
    methods
        function h = MyHandleClass(x)
            h.x = x
        end
        function x = get.x(obj)
        % how is this displayed?
            x = obj.x
        end
    end
    methods (Static)
        function w = my_static_function(z)
        % A static function in :class:`MyHandleClass`.
        %
        % :param z: input z
        % :returns: w

            w = z
        end
    end
end
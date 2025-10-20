classdef FASTLinearizationFile
    %FASTLINEARIZATIONFILE 读取 OpenFAST 线性化输出的 .lin 文件。
    %   该类在 MATLAB 中提供与 Python 版 FASTLinearizationFile 类似的接口。
    %   构造函数可直接给定线性化文件路径，或稍后调用 READ 方法装载数据。
    %
    %   主要属性：
    %     A, B, C, D   - 状态空间矩阵
    %     x, u, y      - 操作点向量
    %     x_info, ...  - 包含坐标系、导数阶数和描述的结构体
    %
    %   主要方法：
    %     read(filename)   - 重新读取指定的 .lin 文件
    %     udescr(), ydescr(), xdescr(), xdotdescr() - 返回简写描述
    %
    %   示例：
    %     lin = FASTLinearizationFile('MyTurbine.1.lin');
    %     A = lin.A; B = lin.B; C = lin.C; D = lin.D;
    %     uDesc = lin.udescr();
    %     yDesc = lin.ydescr();
    %
    %   需要保证包含本文件的目录已加入 MATLAB 路径。

    properties
        filename
        header
        A
        B
        C
        D
        dUdu
        dUdy
        x
        x_info
        xdot
        xdot_info
        u
        u_info
        y
        y_info
        Azimuth
        RotSpeed
        WindSpeed
        JacobiansIncluded logical = false
        NumContStates double = 0
        NumDiscStates double = 0
        NumConstrStates double = 0
        NumInputs double = 0
        NumOutputs double = 0
    end

    methods
        function obj = FASTLinearizationFile(filename)
            if nargin > 0 && ~isempty(filename)
                obj = obj.read(filename);
            end
        end

        function obj = read(obj, filename)
            if nargin > 1 && ~isempty(filename)
                obj.filename = filename;
            end

            if isempty(obj.filename)
                error('FASTLinearizationFile:NoFile', '未指定线性化文件名。');
            end

            fid = fopen(obj.filename, 'r');
            if fid < 0
                error('FASTLinearizationFile:OpenFailed', '无法打开文件: %s', obj.filename);
            end
            cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

            [obj.header, counts, meta] = obj.readHeader(fid);
            obj.NumContStates = counts.nx;
            obj.NumDiscStates = counts.nxd;
            obj.NumConstrStates = counts.nz;
            obj.NumInputs = counts.nu;
            obj.NumOutputs = counts.ny;
            obj.JacobiansIncluded = meta.jacobians;
            obj.Azimuth = meta.azimuth;
            obj.RotSpeed = meta.rotorSpeed;
            obj.WindSpeed = meta.windSpeed;

            % 初始化属性
            obj.A = [];
            obj.B = [];
            obj.C = [];
            obj.D = [];
            obj.dUdu = [];
            obj.dUdy = [];
            obj.x = [];
            obj.x_info = struct('RotatingFrame', {{}}, 'DerivativeOrder', [], 'Description', {{}});
            obj.xdot = [];
            obj.xdot_info = struct('RotatingFrame', {{}}, 'DerivativeOrder', [], 'Description', {{}});
            obj.u = [];
            obj.u_info = struct('RotatingFrame', {{}}, 'DerivativeOrder', [], 'Description', {{}});
            obj.y = [];
            obj.y_info = struct('RotatingFrame', {{}}, 'DerivativeOrder', [], 'Description', {{}});

            while true
                line = fgetl(fid);
                if ~ischar(line)
                    break;
                end
                trimmed = strtrim(line);
                if isempty(trimmed)
                    continue;
                end

                if contains(line, 'Order of continuous states:')
                    [obj.x, obj.x_info] = obj.readOperatingPoints(fid, counts.nx);
                elseif contains(line, 'Order of continuous state derivatives:')
                    [obj.xdot, obj.xdot_info] = obj.readOperatingPoints(fid, counts.nx);
                elseif contains(line, 'Order of inputs')
                    [obj.u, obj.u_info] = obj.readOperatingPoints(fid, counts.nu);
                elseif contains(line, 'Order of outputs')
                    [obj.y, obj.y_info] = obj.readOperatingPoints(fid, counts.ny);
                elseif startsWith(trimmed, 'A:')
                    obj.A = obj.readMatrix(fid, counts.nx, counts.nx, 'A');
                elseif startsWith(trimmed, 'B:')
                    obj.B = obj.readMatrix(fid, counts.nx, counts.nu, 'B');
                elseif startsWith(trimmed, 'C:')
                    obj.C = obj.readMatrix(fid, counts.ny, counts.nx, 'C');
                elseif startsWith(trimmed, 'D:')
                    obj.D = obj.readMatrix(fid, counts.ny, counts.nu, 'D');
                elseif startsWith(trimmed, 'dUdu:')
                    obj.dUdu = obj.readMatrix(fid, counts.nu, counts.nu, 'dUdu');
                elseif startsWith(trimmed, 'dUdy:')
                    obj.dUdy = obj.readMatrix(fid, counts.nu, counts.ny, 'dUdy');
                end
            end
        end

        function names = udescr(obj)
            if ~isempty(obj.u_info) && isfield(obj.u_info, 'Description')
                names = obj.short_descr(obj.u_info.Description);
            else
                names = {};
            end
        end

        function names = ydescr(obj)
            if ~isempty(obj.y_info) && isfield(obj.y_info, 'Description')
                names = obj.short_descr(obj.y_info.Description);
            else
                names = {};
            end
        end

        function names = xdescr(obj)
            if ~isempty(obj.x_info) && isfield(obj.x_info, 'Description')
                names = obj.short_descr(obj.x_info.Description);
            else
                names = {};
            end
        end

        function names = xdotdescr(obj)
            if ~isempty(obj.xdot_info) && isfield(obj.xdot_info, 'Description')
                names = obj.short_descr(obj.xdot_info.Description);
            else
                names = {};
            end
        end
    end

    methods (Access = private)
        function [header, counts, meta] = readHeader(~, fid)
            header = {};
            counts = struct('nx', 0, 'nxd', 0, 'nz', 0, 'nu', 0, 'ny', 0);
            meta = struct('jacobians', false, 'azimuth', NaN, 'rotorSpeed', NaN, 'windSpeed', NaN);

            maxHeaderLines = 200;
            targetFound = false;
            for i = 1:maxHeaderLines
                line = fgetl(fid);
                if ~ischar(line)
                    break;
                end
                header{end+1} = line; %#ok<AGROW>
                if contains(line, 'Jacobians included')
                    targetFound = true;
                    break;
                end
            end

            if ~targetFound
                error('FASTLinearizationFile:FormatError', '未在文件头部找到 ''Jacobians included'' 关键行。');
            end

            counts.nx  = FASTLinearizationFile.safeStr2double(FASTLinearizationFile.extractValue(header, 'Number of continuous states:'), 0);
            counts.nxd = FASTLinearizationFile.safeStr2double(FASTLinearizationFile.extractValue(header, 'Number of discrete states:'), 0);
            counts.nz  = FASTLinearizationFile.safeStr2double(FASTLinearizationFile.extractValue(header, 'Number of constraint states:'), 0);
            counts.nu  = FASTLinearizationFile.safeStr2double(FASTLinearizationFile.extractValue(header, 'Number of inputs:'), 0);
            counts.ny  = FASTLinearizationFile.safeStr2double(FASTLinearizationFile.extractValue(header, 'Number of outputs:'), 0);

            jacStr = FASTLinearizationFile.extractValue(header, 'Jacobians included in this file?');
            meta.jacobians = any(strcmpi(jacStr, {'true', 't', 'yes'}));

            azStr = FASTLinearizationFile.extractValue(header, 'Azimuth:');
            meta.azimuth = FASTLinearizationFile.safeStr2double(azStr, NaN);
            rotStr = FASTLinearizationFile.extractValue(header, 'Rotor Speed:');
            meta.rotorSpeed = FASTLinearizationFile.safeStr2double(rotStr, NaN);
            windStr = FASTLinearizationFile.extractValue(header, 'Wind Speed:');
            meta.windSpeed = FASTLinearizationFile.safeStr2double(windStr, NaN);
        end

        function [op, info] = readOperatingPoints(~, fid, n)
            info = struct('RotatingFrame', {{}}, 'DerivativeOrder', [], 'Description', {{}});
            op = zeros(n, 1);
            if n < 0
                error('FASTLinearizationFile:FormatError', '操作点数量不能为负。');
            end

            % 跳过表头和分隔线
            headerLine = fgetl(fid); %#ok<NASGU>
            separatorLine = fgetl(fid); %#ok<NASGU>

            hasDerivative = contains(headerLine, 'Derivative Order');

            if n == 0
                info.RotatingFrame = {};
                info.DerivativeOrder = [];
                info.Description = {};
                return;
            end

            rotList = strings(n, 1);
            derivList = -ones(n, 1);
            descList = strings(n, 1);

            i = 0;
            while i < n
                line = fgetl(fid);
                if ~ischar(line)
                    error('FASTLinearizationFile:FormatError', '读取操作点时遇到意外的文件结束。');
                end

                if isempty(strtrim(line))
                    continue;
                end

                parts = strsplit(strtrim(line));
                if numel(parts) < 2
                    error('FASTLinearizationFile:FormatError', '操作点行格式不正确: %s', line);
                end

                valueToken = parts{2};
                if contains(valueToken, ',')
                    valueToken = erase(valueToken, ',');
                    rotIdx = 5;
                else
                    rotIdx = 3;
                end

                currentIndex = i + 1;
                op(currentIndex) = str2double(valueToken);
                if isnan(op(currentIndex))
                    error('FASTLinearizationFile:FormatError', '无法解析操作点数值: %s', valueToken);
                end

                if numel(parts) < rotIdx
                    error('FASTLinearizationFile:FormatError', '缺少旋转坐标信息: %s', line);
                end

                rotList(currentIndex) = parts{rotIdx};

                if hasDerivative
                    if numel(parts) < rotIdx + 1
                        error('FASTLinearizationFile:FormatError', '缺少导数阶数: %s', line);
                    end
                    derivList(currentIndex) = str2double(parts{rotIdx + 1});
                    descStart = rotIdx + 2;
                else
                    descStart = rotIdx + 1;
                end

                if descStart <= numel(parts)
                    descList(currentIndex) = strjoin(parts(descStart:end), ' ');
                else
                    descList(currentIndex) = "";
                end

                i = currentIndex;
            end

            info.RotatingFrame = cellstr(rotList);
            info.DerivativeOrder = double(derivList);
            info.Description = cellstr(descList);
        end

        function mat = readMatrix(~, fid, nRows, nCols, name)
            if nRows == 0 || nCols == 0
                mat = zeros(nRows, nCols);
                return;
            end

            mat = zeros(nRows, nCols);
            for r = 1:nRows
                line = fgetl(fid);
                while ischar(line) && isempty(strtrim(line))
                    line = fgetl(fid);
                end
                if ~ischar(line)
                    error('FASTLinearizationFile:FormatError', '矩阵 %s 数据不完整。', name);
                end

                values = sscanf(line, '%f');
                while numel(values) < nCols
                    nextLine = fgetl(fid);
                    if ~ischar(nextLine)
                        error('FASTLinearizationFile:FormatError', '矩阵 %s 的行长度不足。', name);
                    end
                    extra = sscanf(nextLine, '%f');
                    values = [values; extra]; %#ok<AGROW>
                end

                mat(r, :) = values(1:nCols).';
            end

            if any(isnan(mat(:)))
                error('FASTLinearizationFile:FormatError', '矩阵 %s 中存在无法解析的数值。', name);
            end
            if any(isinf(mat(:)))
                error('FASTLinearizationFile:FormatError', '矩阵 %s 中存在无穷大数值。', name);
            end
        end

        function names = short_descr(~, slist)
            if isempty(slist)
                names = {};
                return;
            end

            slist = cellfun(@char, slist, 'UniformOutput', false);
            names = cell(size(slist));
            for i = 1:numel(slist)
                s = string(strtrim(slist{i}));
                if strlength(s) == 0
                    names{i} = '';
                    continue;
                end

                replacements = [
                    "(m/s)", "_[m/s]";
                    "(kW)", "_[kW]";
                    "(deg)", "_[deg]";
                    "(N)", "_[N]";
                    "(kN-m)", "_[kNm]";
                    "(N-m)", "_[Nm]";
                    "(kN)", "_[kN]";
                    "(rpm)", "_[rpm]";
                    "(rad)", "_[rad]";
                    "(rad/s)", "_[rad/s]";
                    "(rad/s^2)", "_[rad/s^2]";
                    "(m/s^2)", "_[m/s^2]";
                    "(deg/s^2)", "_[deg/s^2]";
                    "(m)", "_[m]";
                    ", m/s/s", "_[m/s^2]";
                    ", m/s^2", "_[m/s^2]";
                    ", m/s", "_[m/s]";
                    ", m", "_[m]";
                    ", rad/s/s", "_[rad/s^2]";
                    ", rad/s^2", "_[rad/s^2]";
                    ", rad/s", "_[rad/s]";
                    ", rad", "_[rad]";
                    ", -", "_[-]";
                    ", Nm/m", "_[Nm/m]";
                    ", Nm", "_[Nm]";
                    ", N/m", "_[N/m]";
                    ", N", "_[N]";
                    "(1)", "1";
                    "(2)", "2";
                    "(3)", "3"
                ];

                for k = 1:size(replacements, 1)
                    s = replace(s, replacements(k, 1), replacements(k, 2));
                end

                s = regexprep(s, '\\([^)]*\\)', '');
                s = replace(s, 'ED ', '');
                s = replace(s, 'BD_', 'BD_B');
                s = replace(s, 'IfW ', '');
                s = replace(s, 'Extended input: ', '');
                s = replace(s, '1st tower ', 'qt1');
                s = replace(s, '2nd tower ', 'qt2');

                nd = count(s, 'First time derivative of ');
                if nd > 0
                    s = replace(s, 'First time derivative of ', '');
                    if nd == 1
                        s = "d_" + strtrim(s);
                    elseif nd >= 2
                        s = "dd_" + strtrim(s);
                    end
                end

                s = replace(s, 'Variable speed generator DOF ', 'psi_rot');
                s = replace(s, 'fore-aft bending mode DOF ', 'FA');
                s = replace(s, 'side-to-side bending mode DOF', 'SS');
                s = replace(s, 'bending-mode DOF of blade ', '');
                s = replace(s, ' rotational-flexibility DOF, rad', '-ROT');
                s = replace(s, 'rotational displacement in ', 'rot');
                s = replace(s, 'Drivetrain', 'DT');
                s = replace(s, 'translational displacement in ', 'trans');
                s = replace(s, 'finite element node ', 'N');
                s = replace(s, '-component position of node ', 'posN');
                s = replace(s, '-component inflow on tower node', 'TwrN');
                s = replace(s, '-component inflow on blade 1, node', 'Bld1N');
                s = replace(s, '-component inflow on blade 2, node', 'Bld2N');
                s = replace(s, '-component inflow on blade 3, node', 'Bld3N');
                s = replace(s, '-component inflow velocity at node', 'N');
                s = replace(s, 'X translation displacement, node', 'TxN');
                s = replace(s, 'Y translation displacement, node', 'TyN');
                s = replace(s, 'Z translation displacement, node', 'TzN');
                s = replace(s, 'X translation velocity, node', 'TVxN');
                s = replace(s, 'Y translation velocity, node', 'TVyN');
                s = replace(s, 'Z translation velocity, node', 'TVzN');
                s = replace(s, 'X translation acceleration, node', 'TAxN');
                s = replace(s, 'Y translation acceleration, node', 'TAyN');
                s = replace(s, 'Z translation acceleration, node', 'TAzN');
                s = replace(s, 'X orientation angle, node', 'RxN');
                s = replace(s, 'Y orientation angle, node', 'RyN');
                s = replace(s, 'Z orientation angle, node', 'RzN');
                s = replace(s, 'X rotation velocity, node', 'RVxN');
                s = replace(s, 'Y rotation velocity, node', 'RVyN');
                s = replace(s, 'Z rotation velocity, node', 'RVzN');
                s = replace(s, 'X rotation acceleration, node', 'RAxN');
                s = replace(s, 'Y rotation acceleration, node', 'RAyN');
                s = replace(s, 'Z rotation acceleration, node', 'RAzN');
                s = replace(s, 'X force, node', 'FxN');
                s = replace(s, 'Y force, node', 'FyN');
                s = replace(s, 'Z force, node', 'FzN');
                s = replace(s, 'X moment, node', 'MxN');
                s = replace(s, 'Y moment, node', 'MyN');
                s = replace(s, 'Z moment, node', 'MzN');
                s = replace(s, 'FX', 'Fx');
                s = replace(s, 'FY', 'Fy');
                s = replace(s, 'FZ', 'Fz');
                s = replace(s, 'MX', 'Mx');
                s = replace(s, 'MY', 'My');
                s = replace(s, 'MZ', 'Mz');
                s = replace(s, 'FKX', 'FKx');
                s = replace(s, 'FKY', 'FKy');
                s = replace(s, 'FKZ', 'FKz');
                s = replace(s, 'MKX', 'MKx');
                s = replace(s, 'MKY', 'MKy');
                s = replace(s, 'MKZ', 'MKz');
                s = replace(s, 'Nodes motion', '');
                s = replace(s, 'cosine', 'cos');
                s = replace(s, 'sine', 'sin');
                s = replace(s, 'collective', 'coll.');
                s = replace(s, 'Blade', 'B');
                s = replace(s, 'rotZ', 'TORS-R');
                s = replace(s, 'transX', 'FLAP-D');
                s = replace(s, 'transY', 'EDGE-D');
                s = replace(s, 'rotX', 'EDGE-R');
                s = replace(s, 'rotY', 'FLAP-R');
                s = replace(s, 'flapwise', 'FLAP');
                s = replace(s, 'edgewise', 'EDGE');
                s = replace(s, 'horizontal surge translation DOF', 'Surge');
                s = replace(s, 'horizontal sway translation DOF', 'Sway');
                s = replace(s, 'vertical heave translation DOF', 'Heave');
                s = replace(s, 'roll tilt rotation DOF', 'Roll');
                s = replace(s, 'pitch tilt rotation DOF', 'Pitch');
                s = replace(s, 'yaw rotation DOF', 'Yaw');
                s = replace(s, 'vertical power-law shear exponent', 'alpha');
                s = replace(s, 'horizontal wind speed ', 'WS');
                s = replace(s, 'propagation direction', 'WD');
                s = replace(s, ' pitch command', 'pitch');
                s = replace(s, 'HSS_', 'HSS');
                s = replace(s, 'Bld', 'B');
                s = replace(s, 'tower', 'Twr');
                s = replace(s, 'Tower', 'Twr');
                s = replace(s, 'Nacelle', 'Nac');
                s = replace(s, 'Platform', 'Ptfm');
                s = replace(s, 'SrvD', 'SvD');
                s = replace(s, 'Generator torque', 'Qgen');
                s = replace(s, 'coll. blade-pitch command', 'PitchColl');
                s = replace(s, 'wave elevation at platform ref point', 'WaveElevRefPoint');
                s = replace(s, '1)', '1');
                s = replace(s, '2)', '2');
                s = replace(s, '3)', '3');
                s = replace(s, ',', '');
                s = replace(s, ' ', '');

                names{i} = char(strtrim(s));
            end
        end
    end

    methods (Static, Access = private)
        function val = extractValue(lines, key)
            val = '';
            for i = 1:numel(lines)
                line = lines{i};
                idx = strfind(line, key);
                if ~isempty(idx)
                    tail = strtrim(line(idx + length(key):end));
                    tokens = strsplit(tail);
                    if ~isempty(tokens)
                        val = tokens{1};
                        return;
                    end
                end
            end
        end

        function val = safeStr2double(strVal, defaultVal)
            if nargin < 2
                defaultVal = NaN;
            end
            if isempty(strVal)
                val = defaultVal;
                return;
            end
            val = str2double(strVal);
            if isnan(val)
                val = defaultVal;
            end
        end
    end
end

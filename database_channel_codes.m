function [ code, numCodedBits, numBits ] = database_channel_codes( parameters )
%DATABASE_CHANNEL_CODES Stores details of predefined channel codes.
%   Database Summary:
%
%   PCCCg23g35
%   PCCCg23g35p
%   LDPCk1n2large


switch parameters.chancode
    
    % <<<<< 1/3-rate PCCC >>>>>
    case 'PCCCg23g35'
        
        % Constraint Length
        code.constraint_length = 5;
        
        % Binary Generators
        code.FFbin = [1 0 0 1 1]';     % octal generator 23
        code.FBbin = [1 1 1 0 1]';     % octal generator 35
        code.FF = [35, 23];
        code.FB = 35;
        
        % Convolutional code state termination
        code.termA = 'punc';
        code.termB = 'term';
        if strcmp(code.termA,'term'), code.LA = code.constraint_length-1; else code.LA = 0; end
        if strcmp(code.termB,'term'), code.LB = code.constraint_length-1; else code.LB = 0; end
        
        % Determine number of information/coded bits
        if isfield(parameters,'numBits')==1
            numCodedBits = 3*parameters.numBits + 2*code.LA + 2*code.LB;
        elseif isfield(parameters,'numCodedBits_target')==1
            parameters.numBits = floor( (parameters.numCodedBits_target - 2*code.LA - 2*code.LB)/3 );
            numCodedBits = 3*parameters.numBits + 2*code.LA + 2*code.LB;
        else
            error('Define either ''numBits'' or ''numCodedBits_target''.');
        end
        
        % Code rate
        code.rate = 1/(3+2*(code.LA+code.LB)/parameters.numBits);
        
        % Internal interleaver
        rng(0,'twister'); code.intlvr_int = randperm(parameters.numBits);
        
        % Puncturing
        code.puncMatrix = [1 1; 1 1];
        code.LPM = size(code.puncMatrix,2);
        code.P = [reshape([ones(1,parameters.numBits); repmat(code.puncMatrix,1,parameters.numBits/code.LPM)],1,[]),...
            repmat(code.puncMatrix(1,:),1,2*code.LA/code.LPM),...
            repmat(code.puncMatrix(2,:),1,2*code.LB/code.LPM)];
    % ---------------------------------------------------------------------
    
    
    % <<<<< 1/2-rate PCCC >>>>>
    case 'PCCCg23g35p'
        
        % Constraint Length
        code.constraint_length = 5;
        
        % Binary Generators
        code.FFbin = [1 0 0 1 1]';     % octal generator 23
        code.FBbin = [1 1 1 0 1]';     % octal generator 35
        code.FF = [35, 23];
        code.FB = 35;
        
        % Convolutional code state termination
        code.termA = 'punc';
        code.termB = 'term';
        if strcmp(code.termA,'term'), code.LA = code.constraint_length-1; else code.LA = 0; end
        if strcmp(code.termB,'term'), code.LB = code.constraint_length-1; else code.LB = 0; end
        
        % Determine number of information/coded bits
        if isfield(parameters,'numBits')==1
            numCodedBits = 2*parameters.numBits + code.LA + code.LB;
        elseif isfield(parameters,'numCodedBits_target')==1
            parameters.numBits = floor( (parameters.numCodedBits_target - code.LA - code.LB)/2 );
            numCodedBits = 2*parameters.numBits + code.LA + code.LB;
        else
            error('Define either ''numBits'' or ''numCodedBits_target''.');
        end
        
        % Code rate
        code.rate = 1/(2+(code.LA+code.LB)/parameters.numBits);
        
        % Internal interleaver
        rng(0,'twister'); code.intlvr_int = randperm(parameters.numBits);
        
        % Puncturing
        code.puncMatrix = [1 0; 0 1];
        code.LPM = size(code.puncMatrix,2);
        code.P = [reshape([ones(1,parameters.numBits); repmat(code.puncMatrix,1,parameters.numBits/code.LPM)],1,[]),...
            repmat(code.puncMatrix(1,:),1,2*code.LA/code.LPM),...
            repmat(code.puncMatrix(2,:),1,2*code.LB/code.LPM)];
    % ---------------------------------------------------------------------
      
    
    
    % <<<<< 1/2-rate LDPC large >>>>>
    case 'LDPCk1n2large'
        code.k_rate = 1;
        code.n_rate = 2;
        code.blocksize = 'large'; code.blocksizeValue = 541;
        [~,~,code.G,code.H] = codec_config_LDPC_DMPL_unified(code.k_rate,code.n_rate,code.blocksize,2);
        code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Determine number of coded bits
        numBits = code.blocksizeValue*code.k_rate;
        numCodedBits = code.blocksizeValue*code.n_rate;
        
        % Code rate
        code.rate = code.k_rate/code.n_rate;
        
        % Transpose G matrix
        code.G = code.G.';
        
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,1,1));
    % ---------------------------------------------------------------------
    
    
    % <<<<< 1/2-rate LDPC large2 >>>>>
    case 'LDPCk1n2large2'
        code.k_rate = 1;
        code.n_rate = 2;
        code.blocksize = 'large'; code.blocksizeValue = 541;
        [~,~,code.G,code.H] = codec_config_LDPC_DMPL_unified(code.k_rate,code.n_rate,code.blocksize,2);
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
%         code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Determine number of coded bits
        numBits = code.blocksizeValue*code.k_rate*2;
        numCodedBits = code.blocksizeValue*code.n_rate*2;
        
        % Code rate
        code.rate = code.k_rate/code.n_rate;
        
        % Transpose G matrix
        code.G = code.G.';
        
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,2,1));
    % ---------------------------------------------------------------------
    
    
    % <<<<< 1/2-rate LDPC large4 >>>>>
    case 'LDPCk1n2large4'
        code.k_rate = 1;
        code.n_rate = 2;
        code.blocksize = 'large'; code.blocksizeValue = 541;
        [~,~,code.G,code.H] = codec_config_LDPC_DMPL_unified(code.k_rate,code.n_rate,code.blocksize,2);
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        
        % Determine number of coded bits
        numBits = code.blocksizeValue*code.k_rate*4;
        numCodedBits = code.blocksizeValue*code.n_rate*4;
        
        % Code rate
        code.rate = code.k_rate/code.n_rate;
        
        % Transpose G matrix
        code.G = code.G.';
                
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,4,1));
    % ---------------------------------------------------------------------
    
    
    % <<<<< 1/2-rate LDPC large8 >>>>>
    case 'LDPCk1n2large8'
        code.k_rate = 1;
        code.n_rate = 2;
        code.blocksize = 'large'; code.blocksizeValue = 541;
        [~,~,code.G,code.H] = codec_config_LDPC_DMPL_unified(code.k_rate,code.n_rate,code.blocksize,2);
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        
        % Determine number of coded bits
        numBits = code.blocksizeValue*code.k_rate*8;
        numCodedBits = code.blocksizeValue*code.n_rate*8;
        
        % Code rate
        code.rate = code.k_rate/code.n_rate;
        
        % Transpose G matrix
        code.G = code.G.';
                
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,8,1));
    % ---------------------------------------------------------------------
    
    
    % <<<<< 1/2-rate LDPC large16 >>>>>
    case 'LDPCk1n2large16'
        code.k_rate = 1;
        code.n_rate = 2;
        code.blocksize = 'large'; code.blocksizeValue = 541;
        [~,~,code.G,code.H] = codec_config_LDPC_DMPL_unified(code.k_rate,code.n_rate,code.blocksize,2);
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
        
        % Determine number of coded bits
        numBits = code.blocksizeValue*code.k_rate*16;
        numCodedBits = code.blocksizeValue*code.n_rate*16;
        
        % Code rate
        code.rate = code.k_rate/code.n_rate;
        
        % Transpose G matrix
        code.G = code.G.';
                
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,16,1));
    % ---------------------------------------------------------------------
    
    
    % <<<<< 1/3-rate LDPC large2 >>>>>
    case 'LDPCk1n3large2'
        code.k_rate = 1;
        code.n_rate = 3;
        code.blocksize = 'large'; code.blocksizeValue = 541;
        [~,~,code.G,code.H] = codec_config_LDPC_DMPL_unified(code.k_rate,code.n_rate,code.blocksize,2);
        code.G = [code.G, zeros(size(code.G)); zeros(size(code.G)), code.G];
        code.H = [code.H, zeros(size(code.H)); zeros(size(code.H)), code.H];
%         code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Determine number of coded bits
        numBits = code.blocksizeValue*code.k_rate*2;
        numCodedBits = code.blocksizeValue*code.n_rate*2;
        
        % Code rate
        code.rate = code.k_rate/code.n_rate;
        
        % Transpose G matrix
        code.G = code.G.';
                
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,2,1));
    % ---------------------------------------------------------------------
    
    
    % <<<<< 1/4-rate LDPC large >>>>>
    case 'LDPCk1n4large'
        code.k_rate = 1;
        code.n_rate = 4;
        code.blocksize = 'large'; code.blocksizeValue = 541;
        [~,~,code.G,code.H] = codec_config_LDPC_DMPL_unified(code.k_rate,code.n_rate,code.blocksize,2);
        code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Determine number of coded bits
        numBits = code.blocksizeValue*code.k_rate;
        numCodedBits = code.blocksizeValue*code.n_rate;
        
        % Code rate
        code.rate = code.k_rate/code.n_rate;
        
        % Transpose G matrix
        code.G = code.G.';
        
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,1,1));
    % ---------------------------------------------------------------------

    
    % <<<<< 1/2-rate LDPC generated set 1 >>>>>
    case 'LDPCk1n2gen1'
        
        % Code
        load ldpc_set_1
        code.H = H;
        code.G = G;
        clear H G
        code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Number of coded bits
        [numParity, numCodedBits] = size(code.H);
        numBits = numCodedBits - numParity;
        
        % Code rate
        code.k_rate = 1;
        code.n_rate = 2;
        code.rate = code.k_rate/code.n_rate;
                
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,1,1));
    % ---------------------------------------------------------------------
    

    % <<<<< 1/2-rate LDPC generated set 2 >>>>>
    case 'LDPCk1n2gen2'
        
        % Code
        load ldpc_set_2
        code.H = H;
        code.G = G;
        clear H G
        code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Number of coded bits
        [numParity, numCodedBits] = size(code.H);
        numBits = numCodedBits - numParity;
        
        % Code rate
        code.k_rate = 1;
        code.n_rate = 2;
        code.rate = code.k_rate/code.n_rate;
                        
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,2,1));
    % ---------------------------------------------------------------------


    
    % <<<<< 1/2-rate LDPC generated set 4 >>>>>
    case 'LDPCk1n2gen4'
        
        % Code
        load ldpc_set_4
        code.H = H;
        code.G = G;
        clear H G
        code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Number of coded bits
        [numParity, numCodedBits] = size(code.H);
        numBits = numCodedBits - numParity;
        
        % Code rate
        code.k_rate = 1;
        code.n_rate = 2;
        code.rate = code.k_rate/code.n_rate;
                        
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,4,1));
    % ---------------------------------------------------------------------

    
    % <<<<< 1/2-rate LDPC generated set 8 >>>>>
    case 'LDPCk1n2gen8'
        
        % Code
        load ldpc_set_8
        code.H = H;
        code.G = G;
        clear H G
        code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Number of coded bits
        [numParity, numCodedBits] = size(code.H);
        numBits = numCodedBits - numParity;
        
        % Code rate
        code.k_rate = 1;
        code.n_rate = 2;
        code.rate = code.k_rate/code.n_rate;
                        
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,8,1));
    % ---------------------------------------------------------------------


    % <<<<< 1/2-rate LDPC generated set 16 >>>>>
    case 'LDPCk1n2gen16'
        
        % Code
        load ldpc_set_16
        code.H = H;
        code.G = G;
        clear H G
        code.hEnc = comm.LDPCEncoder(sparse(code.H));
        
        % Number of coded bits
        [numParity, numCodedBits] = size(code.H);
        numBits = numCodedBits - numParity;
        
        % Code rate
        code.k_rate = 1;
        code.n_rate = 2;
        code.rate = code.k_rate/code.n_rate;
                        
        % Information bit index
        ind = [zeros(code.blocksizeValue*code.n_rate - code.blocksizeValue*code.k_rate,1);...
               ones(code.blocksizeValue*code.k_rate,1)];
        code.bitIndex = find(repmat(ind,16,1));
    % ---------------------------------------------------------------------




    
end

end


#! /bin/env ruby


################################################################################

# updated on 2023-01-11
#   auto-detect input type (marginal or joint)
# updated on 2023-01-09
#   --rtc a folder as an input allowed
# updated on 2022-11-18
#   bugs fixed
# updated on 2022-11-18 for non-independent ASR based on scm


################################################################################

require 'find'
require 'getoptlong'
require 'csv'
require 'parallel'
require 'bio-nwk'
require 'colorize'


################################################################################
class Dir
    def self.mkdirs(path)
        if(!File.directory?(path))
            if(!mkdirs(File.dirname(path)))
                return false;
            end
            mkdir(path)
        end
        return true
    end
end


################################################################################
def mkdir_with_force(outdir, is_force=false, is_tolerate=false)
  if outdir.class != String
    raise "outdir wrong? Exiting ......"
  end

  if ! Dir.exists?(outdir)
    `mkdir -p #{outdir}`
  else
    if is_tolerate
      ;
    elsif is_force
      `rm -rf #{outdir}`
      `mkdir -p #{outdir}`
    else
      raise "The outdir #{outdir} has already existed!"
    end
  end
end


def read_infiles(indir, suffix='', is_all_subfolder=false)
  infiles = Array.new
  if ! is_all_subfolder
    Dir.foreach(indir) do |b|
      next if b =~ /^\./
      if suffix.is_a?(String)
        if suffix != ''
          next if b !~ /#{suffix}$/
        end
      elsif suffix.is_a?(Array)
        next unless suffix.any?{|i| b =~ /#{i}$/ }
      end
      infiles << File.join(indir, b)
    end
  else
    Find.find(indir) do |path|
      next if File.directory?(path)
      next if File.basename(path) =~ /^\./
      infiles << path if suffix.is_a?(String) ? path =~ /\.#{suffix}$/ : suffix.any?{|i| path =~ /#{i}$/ }
    end
  end
  return(infiles)
end


def getFilesBySuffices(indir, suffices)
  files = Array.new
  infiles = read_infiles(indir)
  infiles.each do |infile|
    if suffices.include?(File.extname(infile))
      files << infile
    end
  end
  return(files)
end


def get_file_path(file)
  path = File.symlink?(file) ? File.readlink(file) : file
  return(path)
end

################################################################################

# all credits to Wolfgang Teuber (https://gitlab.com/knugie)

def weighted_rand(weights = {})
  raise 'Probabilities must sum up to 1' unless weights.values.inject(&:+) == 1.0
  raise 'Probabilities must not be negative' unless weights.values.all? { |p| p >= 0 }
  # Do more sanity checks depending on the amount of trust in the software component using this method,
  # e.g. don't allow duplicates, don't allow non-numeric values, etc.

  # Ignore elements with probability 0
  weights = weights.reject { |k, v| v == 0.0 }   # e.g. => {"a"=>0.4, "b"=>0.4, "c"=>0.2}

  # Accumulate probabilities and map them to a value
  u = 0.0
  ranges = weights.map { |v, p| [u += p, v] }   # e.g. => [[0.4, "a"], [0.8, "b"], [1.0, "c"]]

  # Generate a (pseudo-)random floating point number between 0.0(included) and 1.0(excluded)
  u = rand   # e.g. => 0.4651073966724186

  # Find the first value that has an accumulated probability greater than the random number u
  ranges.find { |p, v| p > u }.last   # e.g. => "b"
end

################################################################################

class Bio::Tree
  def de_no_for_tips
    allTips.each do |tip|
      tip.name = tip.name.split(' ')[1,100].join('_')
    end
  end
end


class Subclade
  attr_accessor :name, :num
  def initialize(arr)
    #@name = arr
    @name = arr[0].is_a?(Array) ? arr.map{|i|i.sort} : arr.sort
    #@name = arr[0] =~ /,/ ? arr.map{|i|i.split(',').sort} : arr.sort
    @num = nil
  end
end


class Symbiont < Subclade
  attr_accessor :hosts
  def initialize(arr)
    @name = arr
    @hosts = Array.new
  end
  def is_co_evolve?(data:, is_strict:)
    @host2prob = Hash.new
    @hosts.map{|i| @host2prob[i]=i.prob }
    selected_host = is_strict ? @hosts.sort_by{|i| i.prob }.reverse[0] : weighted_rand(@host2prob)
    if data[selected_host.num-1].to_f >= data[@num-1].to_f
      return(true)
    else
      return(false)
    end
  end
  def is_co_evolve2?(data:, is_strict:)
    @host2prob = Hash.new
    @hosts.map{|i| @host2prob[i]=i.prob }
    selected_host = is_strict ? @hosts.sort_by{|i| i.prob }.reverse[0] : weighted_rand(@host2prob)
    if selected_host.num.zip(@num).all?{|a, b| data[a-1].to_f >= data[b-1].to_f }
      return(true)
    else
      return(false)
    end
  end
end


class Host < Subclade
  attr_accessor :prob
  def prob
    @prob
  end
  def prob=(prob)
    @prob = prob.to_f
  end
end


##############################################
def get_lca_bootstrap(names, name2node, tree)
  nodes = names.map{|i| name2node[i] }
  lca = tree.lowest_common_ancestor(nodes[0], nodes[1])
  begin
    return(lca.bootstrap)
  rescue
    raise nodes.join("\t")
  end
end


def read_mcmctree_out(file)
  is_start = false
  tree = nil
  in_fh = File.open(file, 'r')
  #(((1_t5, 2_t9) 33 , ((3_t21, (4_t22, 5_t18) 36 ) 35 , (6_t3, 7_t6) 37 ) 34 ) 32 , (((((8_t8, 9_t12) 42 , (10_t23, 11_t14) 43 ) 41 , (12_t2, 13_t4) 44 ) 40 , 14_t24) 39 , ((((15_t30, 16_t20) 48 , 17_t13) 47 , (((18_t11, 19_t10) 51 , ((20_t1, 21_t26) 53 , (22_t7, (23_t19, 24_t27) 55 ) 54 ) 52 ) 50 , ((((25_t25, 26_t17) 59 , 27_t15) 58 , 28_t16) 57 , 29_t29) 56 ) 49 ) 46 , 30_t28) 45 ) 38 ) 31 ;
  in_fh.each_line do |line|
    line.chomp!
    if is_start
      tree = getTreeObjFromNwkString(line)
      break
    end
    is_start = true if line =~ /^Species tree for FigTree/
  end
  in_fh.close

  tree.de_no_for_tips

  #get_internal_node_index(tree)
  return(tree)
end


def get_rtc(file, root_two_children_names)
  rtc_info = Hash.new
  in_fh = File.open(file, 'r')
  in_fh.each_line do |line|
    line.chomp!
    next if line =~ /^#|^$/
    line_arr = line.split("\t")

    symbiont = line_arr[0].gsub('_', '_').split(',')
    rtc = Symbiont.new(symbiont)

    (1..line_arr.size-1).each do |index|
      ele = line_arr[index]
      host = ele.split(':')[0].gsub('_', '_').split(',')
      prob = ele.split(':')[1]

      host_obj = Host.new(host)
      host_obj.prob = prob
      rtc.hosts << host_obj
    end

    root_obj = Host.new(root_two_children_names)

    curr_total_prob = rtc.hosts.map{|host|host.prob}.reduce(&:+)
    if curr_total_prob < 1
      root_obj.prob = 1 - curr_total_prob
      rtc.hosts << root_obj
    end

    rtc_info[symbiont] = rtc
  end
  in_fh.close
  return(rtc_info)
end


def get_scm(file, root_two_children_names)
  scm_info = Hash.new
  symbionts = nil
  in_fh = File.open(file, 'r')
  rtc = nil

  in_fh.each_line do |line|
    line.chomp!
    line_arr = line.split("\t")
    next if line =~ /^#|^$/
    if $. == 1
      symbionts = line_arr.map{|i| i.gsub(' ', '_').split(',') }
      rtc = Symbiont.new(symbionts)
      next
    end

    prob = line_arr[-1].to_f
    hosts = line_arr[0, line_arr.size-1].map{|i|i.split(',')}
    hosts.map!{|a| a[0] == 'root' ? root_two_children_names : a } #2023-01, in case of >=1 root (fl; Z)
    host_obj = Host.new(hosts)
    host_obj.prob = prob
    rtc.hosts << host_obj
  end

  root_obj = Host.new([root_two_children_names] * symbionts.size)
  curr_total_prob = rtc.hosts.map{|host|host.prob}.reduce(&:+)
  if curr_total_prob < 1
    root_obj.prob = 1 - curr_total_prob
    rtc.hosts << root_obj
  end
  scm_info[symbionts] = rtc

  in_fh.close
  return(scm_info)
end


##############################################
def auto_detect_mj(rtc_files) # identify whether marginal or joint
  first_file = rtc_files[0]
  in_fh = File.open(first_file, 'r')
  first_line = in_fh.readline.chomp
  if first_line =~ /:[10] (\b | \.\d+)/x
    type = 'marginal'
  else
    type = 'joint'
  end
  in_fh.close
  STDERR.puts "type auto-detected\t" + type.colorize(:yellow)
  return(type)
end


##############################################
if __FILE__ == $0
  indir = nil
  infile = nil
  type = nil
  rtc_files = Array.new
  scm_file = nil
  scm_files = Array.new
  mcmctxt_file = nil
  is_rtc = true
  is_renum = true
  is_strict = false

  rtc_info = Hash.new
  scm_info = Hash.new

  ##############################################
  opts = GetoptLong.new(
    ['--indir', GetoptLong::REQUIRED_ARGUMENT],
    ['-i', GetoptLong::REQUIRED_ARGUMENT],
    ['--marginal', GetoptLong::REQUIRED_ARGUMENT],
    ['--joint', GetoptLong::REQUIRED_ARGUMENT],
    ['--rtc', '--rrtc', GetoptLong::REQUIRED_ARGUMENT],
    ['--mcmctxt', GetoptLong::REQUIRED_ARGUMENT],
    ['--is_rtc', '--is_rrtc', '--rrtc_file', '--rtc_file', GetoptLong::REQUIRED_ARGUMENT],
    ['--strict', '--is_strict', GetoptLong::NO_ARGUMENT],
    ['--no_renum', '--no_re_num', GetoptLong::NO_ARGUMENT],
  )


  opts.each do |opt, value|
    case opt
      when '--indir'
        indir = value
      when '-i'
        infile = value
      when '--marginal'
        rtc_files = File.directory?(value) ? read_infiles(value) : [value]
        type = 'marginal'
      when '--joint'
        rtc_files = File.directory?(value) ? read_infiles(value) : [value]
        type = 'joint'
      when /^--r?rtc(file)?$/
        rtc_files = File.directory?(value) ? read_infiles(value) : [value]
      when '--mcmctxt'
        mcmctxt_file = value
      when '--is_rtc', '--is_rrtc'
        is_rtc = value =~ /^true|T$/i ? true : false
      when '--is_strict', '--strict'
        is_strict = true
        STDERR.puts "is_strict:\t" + "true".colorize(:yellow)
      when '--no_renum', '--no_re_num'
        is_renum = false
    end
  end


  ##############################################
  unless indir.nil?
    infile = File.join(indir, 'out') # file 'out' from the mcmctree output
    mcmctxt_file = File.join(indir, 'mcmc.txt')
    if is_output
      outdir = File.join(File.dirname(mcmctxt_file), '')
      out_fh = File.open(outfile, 'w')
    end
  end

  type = auto_detect_mj(rtc_files) if type.nil? # identify whether marginal or joint


  ##############################################
  tree = read_mcmctree_out(infile) # "out", NOT "mcmc.txt"

  root_two_children_names = tree.children(tree.root).map{|i|tree.tips(i)[0].name}

  case type
    when 'marginal'
      #rtc_info = get_rtc(rtc_file, root_two_children_names)
      rtc_files.each do |rtc_file|
        rtc_info.merge! get_rtc(rtc_file, root_two_children_names)
      end
    when 'joint'
    rtc_files.each do |scm_file|
      # note scm_info
      scm_info.merge! get_scm(scm_file, root_two_children_names)
    end
  else
    raise "rtc_file or scm_file has to be provided! Exiting ......"
  end

  name2node, node2name = tree.getNameNodeRela
  rtc_info.delete_if{|names, rtc| not names.all?{|name|name2node.include?(name)}} unless rtc_info.nil?
  scm_info.delete_if{|names, rtc| not names.flatten.all?{|name|name2node.include?(name)}} unless rtc_info.nil?

  minus = tree.allTips.size - 1


  ##############################################
  rtc_info.each_pair do |names, rtc|
    begin
      rtc.num = get_lca_bootstrap(names, name2node, tree) - minus
      rtc.hosts.map{|obj|obj.num = get_lca_bootstrap(obj.name, name2node, tree) - minus}
    rescue
      raise "species #{names} or #{rtc.hosts} not found"
    end
  end

  scm_info.each_pair do |names, rtc|
    begin
      rtc.num = names.map{|names2|get_lca_bootstrap(names2, name2node, tree) - minus}
      scm_info[names].num = rtc.num
      rtc.hosts.each do |obj|
        obj_nums = Array.new
        names = obj.name
        names.each do |names|
          obj_nums << get_lca_bootstrap(names, name2node, tree) - minus
        end
        obj.num = obj_nums
      end
    end
  end


  ##############################################
  headers = CSV.open(mcmctxt_file, &:readline)

  col_data = Array.new
  CSV.foreach(mcmctxt_file) do |row|
    data = row[0].split("\t")
    if is_rtc
      if not rtc_info.empty?
        #col_data << row if rtc_info.all?{|rtc_name, rtcs| rtcs.all?{|rtc| rtc.is_co_evolve?(data)} }
        col_data << row if rtc_info.all?{|rtc_name, rtc| rtc.is_co_evolve?(data:data, is_strict:is_strict) }
      elsif not scm_info.empty?
        col_data << row if scm_info.all?{|rtc_name, rtc| rtc.is_co_evolve2?(data:data, is_strict:is_strict) }
      end
    else
      col_data << row
    end
  end

  STDERR.puts "# of samples after filtering\t" + col_data.size.to_s.colorize(:red)

  col_data.each_with_index do |row_arr, index| # row_arr: single-element array
    if is_renum and index > 0
      posteriors = row_arr[0].split("\t")
      posteriors[0] = index
      row_arr = [posteriors.join("\t")]
    end
    puts row_arr
  end
end



FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/latch-base:02ab-main

# Its easy to build binaries from source that you can later reference as
# subprocesses within your workflow.
# RUN curl -L https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.4/bowtie2-2.4.4-linux-x86_64.zip/download -o bowtie2-2.4.4.zip &&\
# unzip bowtie2-2.4.4.zip &&\
# mv bowtie2-2.4.4-linux-x86_64 bowtie2

# Or use managed library distributions through the container OS's package
# manager.
# RUN apt-get update -y &&\
#  apt-get install -y autoconf samtools


# Installing utilities
RUN  apt-get update -y &&\
  apt-get install -y wget &&\
  apt-get install -y sudo 

# Installing java

# RUN sudo apt-get install -y default-jre &&\
# apt-get install -y default-jdk
RUN apt-get update \
  && apt-get install -y default-jre-headless

# Install trimmomatric from source
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip &&\
  unzip Trimmomatic-0.39.zip &&\
  cd Trimmomatic-0.39 &&\
  mv trimmomatic-0.39.jar $HOME
COPY trimmomatic $HOME
RUN chmod +x $HOME/trimmomatic



# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
RUN python3 -m pip install --upgrade latch
WORKDIR /root

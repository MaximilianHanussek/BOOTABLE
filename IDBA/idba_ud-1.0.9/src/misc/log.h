/**
 * @file log.h
 * @brief LogThread Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.8
 * @date 2011-12-23
 */

#ifndef __MISC_LOG_H_

#define __MISC_LOG_H_

#include <sys/types.h>
#include <pthread.h>

#include <cstdio>
#include <string>

class LogProcess
{
public:
    LogProcess(const std::string &log_file);
    ~LogProcess();
    
private:
    pid_t pid;
};

class LogThread
{
public:
    LogThread(const std::string &log_file);
    ~LogThread();

private:
    static void *LogThreadFunc(void *p);

    int in_fd;
    int out_fd;
    std::string log_file_;
    pthread_t thread_;
};

#endif


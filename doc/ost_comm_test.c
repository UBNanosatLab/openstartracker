#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>

#define BUFFER_SIZE 1024
#define MAX_STARS 50

/* Coordinates of a star in both camera and ECI frames, and its weight */
typedef struct {
    double cam_x;
    double cam_y;
    double cam_z;
    double eci_x;
    double eci_y;
    double eci_z;
    double weight;
} star_t;

/* sendImageCommand:
 *   This function establishes a TCP connection to a server, sends an image command,
 *   and receives a response from the server.
 *
 * Input arguments:
 *  - server_ip: IP address of the server.
 *  - server_port: Port number of the server.
 *  - command: Image command to send to the server.
 *  - output_buffer: Starer to the output buffer where the response will be stored.
 *  - max_output_size: Maximum size of the output buffer to prevent buffer overflow.
 *
 * Output arguments:
 *  - output_buffer: The response from the server will be stored in this buffer.
 *
 * Return value:
 *  - If the function succeeds, it returns the number of bytes received and stored in the output buffer.
 *  - If an error occurs, it returns -1.
 */
int sendImageCommand(const char* server_ip, int server_port, const char* command, char* output_buffer, int max_output_size) {
    int sockfd;
    struct sockaddr_in server_addr;

    // Create socket
    if ((sockfd = socket(AF_INET, SOCK_STREAM, 0)) == -1) {
        perror("Socket creation failed");
        return -1;
    }

    // Fill server information
    server_addr.sin_family = AF_INET;
    server_addr.sin_port = htons(server_port);
    server_addr.sin_addr.s_addr = inet_addr(server_ip);

    // Connect to the server
    if (connect(sockfd, (struct sockaddr *)&server_addr, sizeof(server_addr)) != 0) {
        perror("Connection to server failed");
        close(sockfd);
        return -1;
    }

    // Send the image command to the server
    if (send(sockfd, command, strlen(command), 0) == -1) {
        perror("Send failed");
        close(sockfd);
        return -1;
    }

    // Receive response from the server
    int bytes_received = recv(sockfd, output_buffer, max_output_size - 1, 0);
    if (bytes_received == -1) {
        perror("Receive failed");
        close(sockfd);
        return -1;
    }
    output_buffer[bytes_received] = '\0'; // Null-terminate the received data

    // Close the socket
    close(sockfd);

    // Return the number of bytes received
    return bytes_received;
}

/* parse_ost_string:
 *   Parses a openstartracker string into an array of star structures.
 *
 * Input arguments:
 *  - ost_str: The openstartracker string to parse.
 *      Example (2 stars): "1.0,2.0,3.0,4.0,5.0,6.0,7.0 2.0,3.0,4.0,5.0,6.0,7.0,8.0"
 *  - max_num_stars: The maximum number of stars to parse.
 *      If the openstartracker string contains more stars than max_num_stars, only the first max_num_stars will be parsed.
 * 
 * Output arguments:
 *  - stars: The parsed stars will be stored in this array.
 *      The array must be pre-allocated by the caller.
 *      The array must have at least max_num_stars elements.
 *      If the openstartracker string contains fewer stars than max_num_stars, the remaining elements will be left unchanged.
 *
 * Return value:
 *  - The number of stars successfully parsed from the openstartracker string.
 */
int parse_ost_string(const char *ost_str, int max_num_stars, star_t *stars) {
    char *ptr = (char *)ost_str;
    int i = 0;
    for (; i < max_num_stars; i++) { 
        // Parse the openstartracker string into a star structure
        int n = sscanf(ptr, "%lf,%lf,%lf,%lf,%lf,%lf,%lf", &stars[i].cam_x, &stars[i].cam_y, &stars[i].cam_z, &stars[i].eci_x, &stars[i].eci_y, &stars[i].eci_z, &stars[i].weight);
        if (n != 7) {
            break;
        }
        // Move the pointer to the next star
        ptr = strchr(ptr, ' ');
        if (ptr == NULL) {
            break;
        }
        // Skip the space character
        ptr++;
    }
    return i;
}

star_t stars[MAX_STARS];
int main() {
    const char *server_ip = "127.0.0.1";  // IP address of the server
    int server_port = 8010;               // Port number the server is listening on
    int max_num_stars = MAX_STARS;       // Maximum number of stars to parse

    // Send the image command to the server
    char output_buffer[BUFFER_SIZE];
    int bytes_received = sendImageCommand(server_ip, server_port, "rgb.solve_image('science_cam_may8_0.05sec_gain40/samples/img0.png')", output_buffer, BUFFER_SIZE);
    if (bytes_received == -1) {
        return 1;
    }

    // Parse the openstartracker string into an array of star structures
    int num_stars = parse_ost_string(output_buffer, max_num_stars, stars);
    if (num_stars == 0) {
        printf("No stars found\n");
        return 1;
    }

    // Print the stars
    for (int i = 0; i < num_stars; i++) {
        printf("Star %d: cam=(%f, %f, %f), eci=(%f, %f, %f), weight=%f\n", i, stars[i].cam_x, stars[i].cam_y, stars[i].cam_z, stars[i].eci_x, stars[i].eci_y, stars[i].eci_z, stars[i].weight);
    }

    return 0;
}
